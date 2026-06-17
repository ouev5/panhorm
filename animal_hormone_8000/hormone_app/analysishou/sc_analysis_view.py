# ！！！关键：在所有导入前，通过Numba API强制禁用缓存（比配置项优先级更高）
import numba
numba.config.DISABLE_CACHE = True  # 直接禁用JIT缓存，比config更底层

# 解决Matplotlib权限问题
import os
mpl_cache_dir = "/tmp/matplotlib_cache"
os.makedirs(mpl_cache_dir, exist_ok=True)
os.environ["MPLCONFIGDIR"] = mpl_cache_dir

# 强制Matplotlib后端
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 其他导入（保持不变）
import uuid
import tempfile
import zipfile
import logging
import threading
from datetime import datetime
from django.http import JsonResponse, HttpResponse
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from django.conf import settings
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from django.db import models
from django.contrib.auth.models import User
import pandas as pd
import numpy as np

# ！！！延迟导入scanpy，确保Numba配置完全生效
import scanpy as sc

# 设置日志
logger = logging.getLogger(__name__)

# 媒体目录设置（保持不变）
sc_uploads_dir = os.path.join(settings.MEDIA_ROOT, 'sc_uploads')
os.makedirs(sc_uploads_dir, exist_ok=True, mode=0o755)
sc_results_dir = os.path.join(settings.MEDIA_ROOT, 'sc_results')
os.makedirs(sc_results_dir, exist_ok=True, mode=0o755)


# 数据模型（保持不变）
class ScAnalysisTask(models.Model):
    TASK_STATUS = (
        ('pending', '等待中'),
        ('running', '运行中'),
        ('completed', '已完成'),
        ('failed', '失败')
    )
    
    task_id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    user = models.ForeignKey(User, on_delete=models.CASCADE, null=True, blank=True)
    task_name = models.CharField(max_length=200, default="未命名分析")
    status = models.CharField(max_length=20, choices=TASK_STATUS, default='pending')
    progress = models.IntegerField(default=0)
    created_at = models.DateTimeField(auto_now_add=True)
    completed_at = models.DateTimeField(null=True, blank=True)
    error_message = models.TextField(null=True, blank=True)
    
    qc_params = models.JSONField(default=dict, blank=True)
    normalization_params = models.JSONField(default=dict, blank=True)
    clustering_params = models.JSONField(default=dict, blank=True)
    annotation_params = models.JSONField(default=dict, blank=True)
    de_params = models.JSONField(default=dict, blank=True)
    enrichment_params = models.JSONField(default=dict, blank=True)
    
    def __str__(self):
        return f"{self.task_name} ({self.task_id})"


class ScUploadedFile(models.Model):
    FILE_TYPES = (
        ('sc-data', '单细胞表达数据'),
        ('metadata', '细胞元数据'),
        ('gene-list', '基因列表'),
        ('sample-info', '样本信息'),
        ('marker-genes', '标记基因'),
        ('custom-genesets', '自定义基因集')
    )
    
    file_id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    task = models.ForeignKey(ScAnalysisTask, related_name='uploaded_files', on_delete=models.CASCADE)
    file_type = models.CharField(max_length=20, choices=FILE_TYPES)
    file = models.FileField(upload_to='sc_uploads/')
    filename = models.CharField(max_length=255)
    file_size = models.IntegerField(help_text="文件大小(字节)")
    uploaded_at = models.DateTimeField(auto_now_add=True)
    
    def __str__(self):
        return f"{self.get_file_type_display()}: {self.filename}"


class ScAnalysisResult(models.Model):
    RESULT_TYPES = (
        ('anndata', '处理后的AnnData文件'),
        ('de-genes', '差异表达基因'),
        ('cell-annotations', '细胞注释'),
        ('enrichment', '富集分析结果'),
        ('visualizations', '可视化图表'),
        ('report', '分析报告')
    )
    
    result_id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    task = models.ForeignKey(ScAnalysisTask, related_name='results', on_delete=models.CASCADE)
    result_type = models.CharField(max_length=20, choices=RESULT_TYPES)
    file = models.FileField(upload_to='sc_results/')
    filename = models.CharField(max_length=255)
    created_at = models.DateTimeField(auto_now_add=True)
    
    def __str__(self):
        return f"{self.get_result_type_display()}: {self.filename}"


# 分析工具类（增强版）
class ScAnalysisTool:
    def __init__(self, task):
        self.task = task
        self.adata = None
        self.results_dir = os.path.join(sc_results_dir, str(task.task_id))
        os.makedirs(self.results_dir, exist_ok=True, mode=0o755)
    
    def update_progress(self, progress, status_text=None):
        self.task.progress = progress
        if status_text:
            self.task.status = status_text
        self.task.save()
        logger.info(f"Task {self.task.task_id} progress updated to {progress}%")
    
    def load_data(self, sc_data_file):
        self.update_progress(5, 'running')
        logger.info(f"Loading data from {sc_data_file}")
        
        file_ext = os.path.splitext(sc_data_file)[1].lower()
        if file_ext == '.h5ad':
            self.adata = sc.read_h5ad(sc_data_file)
        elif file_ext == '.csv':
            df = pd.read_csv(sc_data_file, index_col=0)
            self.adata = sc.AnnData(df.T)
        elif file_ext == '.mtx':
            self.adata = sc.read_mtx(sc_data_file).T
        elif file_ext == '.loom':
            self.adata = sc.read_loom(sc_data_file)
        else:
            raise ValueError(f"不支持的文件格式: {file_ext}，支持格式：.h5ad, .csv, .mtx, .loom")
        
        logger.info(f"Data loaded successfully. Shape: {self.adata.shape}")
        return self.adata
    
    def quality_control(self, params):
        self.update_progress(10, 'running')
        logger.info("Starting quality control...")
        logger.info(f"Quality Control Params Received: {params}") # 添加日志
        
        # 计算线粒体基因比例
        # 尝试匹配线粒体基因名称，不同物种的前缀可能不同
        mt_prefixes = ['MT-', 'mt-', 'Mt-'] # 可以根据需要添加更多前缀
        mt_gene_mask = pd.Series(self.adata.var_names).str.startswith(tuple(mt_prefixes))
        self.adata.var['mt'] = mt_gene_mask.values
        sc.pp.calculate_qc_metrics(
            self.adata, 
            qc_vars=['mt'], 
            percent_top=None, 
            log1p=False, 
            inplace=True
        )
        
        # 过滤细胞 - 只使用一个条件
        min_genes_per_cell = params.get('minGenesPerCell')
        max_genes_per_cell = params.get('maxGenesPerCell')
        
        # --- 添加调试日志 ---
        logger.info(f"Applying cell filter - min_genes: {min_genes_per_cell}, max_genes: {max_genes_per_cell}")
        # --- 添加调试日志 ---
        
        # 确保 min 和 max 不会同时被设置
        if min_genes_per_cell is not None and max_genes_per_cell is not None:
            logger.warning("Warning: Both minGenesPerCell and maxGenesPerCell are set. Only applying minGenesPerCell.")
            # 或者你可以选择报错
            # raise ValueError("Cannot set both minGenesPerCell and maxGenesPerCell at the same time.")
        
        if min_genes_per_cell is not None:
            logger.info(f"Filtering cells with min_genes >= {min_genes_per_cell}")
            sc.pp.filter_cells(self.adata, min_genes=min_genes_per_cell)
        elif max_genes_per_cell is not None:
            logger.info(f"Filtering cells with max_genes <= {max_genes_per_cell}")
            sc.pp.filter_cells(self.adata, max_genes=max_genes_per_cell)
        # 如果两者都为 None，则不进行基于基因数量的细胞过滤
        
        # 过滤基因 - 只使用一个条件
        min_cells_per_gene = params.get('minCellsPerGene')
        max_counts_per_gene = params.get('maxCountsPerGene') # 假设你也有这个参数
        
        # --- 添加调试日志 ---
        logger.info(f"Applying gene filter - min_cells: {min_cells_per_gene}, max_counts: {max_counts_per_gene}")
        # --- 添加调试日志 ---
        
        if min_cells_per_gene is not None and max_counts_per_gene is not None:
            logger.warning("Warning: Both minCellsPerGene and maxCountsPerGene are set. Only applying minCellsPerGene.")
            # 或者你可以选择报错
            # raise ValueError("Cannot set both minCellsPerGene and maxCountsPerGene at the same time.")

        if min_cells_per_gene is not None:
            logger.info(f"Filtering genes with min_cells >= {min_cells_per_gene}")
            sc.pp.filter_genes(self.adata, min_cells=min_cells_per_gene)
        elif max_counts_per_gene is not None:
            logger.info(f"Filtering genes with max_counts <= {max_counts_per_gene}")
            sc.pp.filter_genes(self.adata, max_counts=max_counts_per_gene)
        # 如果两者都为 None，则不进行基于细胞数量的基因过滤

        # 过滤高线粒体含量的细胞
        max_mito_ratio = params.get('maxMitoRatio', 10) / 100 # 默认 10%
        logger.info(f"Filtering cells with mitochondrial ratio < {max_mito_ratio * 100}%")
        self.adata = self.adata[self.adata.obs.pct_counts_mt < max_mito_ratio, :]
        
        logger.info(f"After QC: {self.adata.shape[0]} cells, {self.adata.shape[1]} genes")
        self.update_progress(20, 'running')
        return self.adata
    
    def normalization(self, params):
        self.update_progress(25, 'running')
        logger.info("Starting normalization...")
        
        norm_method = params.get('normMethod', 'log1p')
        
        if norm_method == 'log1p':
            logger.info("Applying log1p normalization...")
            sc.pp.normalize_total(self.adata, target_sum=1e4)
            sc.pp.log1p(self.adata)
        elif norm_method == 'sctransform':
            try:
                import sctransform
                logger.info("Applying SCTransform normalization...")
                model = sctransform.SCTransform(self.adata.X.T)
                self.adata.X = model.transform().T
            except ImportError:
                raise ImportError("使用sctransform需要先安装：pip install sctransform")
        elif norm_method == 'rle':
            logger.info("Applying RLE normalization...")
            sc.pp.normalize_total(self.adata, target_sum=1e4)
        
        # 选择高变基因
        n_hvg = params.get('nHVG', 2000)
        logger.info(f"Selecting top {n_hvg} highly variable genes...")
        sc.pp.highly_variable_genes(
            self.adata, 
            min_mean=0.0125, 
            max_mean=3, 
            min_disp=0.5, 
            n_top_genes=min(n_hvg, self.adata.shape[1])  # 确保不超过基因总数
        )
        self.adata = self.adata[:, self.adata.var.highly_variable]
        
        logger.info(f"After normalization: {self.adata.shape[0]} cells, {self.adata.shape[1]} HVGs")
        self.update_progress(35, 'running')
        return self.adata
    
    def dimensionality_reduction(self, params):
        self.update_progress(40, 'running')
        logger.info("Starting dimensionality reduction...")
        
        logger.info("Scaling data...")
        sc.pp.scale(self.adata, max_value=10)
        
        n_pca = min(params.get('nPcaComponents', 50), self.adata.shape[1]-1)  # 确保不超过特征数
        logger.info(f"Computing PCA with {n_pca} components...")
        sc.tl.pca(self.adata, svd_solver='arpack', n_comps=n_pca)
        
        n_neighbors = min(params.get('umapNeighbors', 15), self.adata.shape[0]-1)
        logger.info(f"Computing neighbors with n_neighbors={n_neighbors}...")
        sc.pp.neighbors(
            self.adata, 
            n_neighbors=n_neighbors, 
            n_pcs=n_pca
        )
        
        logger.info("Computing UMAP...")
        sc.tl.umap(self.adata, min_dist=params.get('umapMinDist', 0.5))
        
        cluster_method = params.get('clusterMethod', 'leiden')
        resolution = params.get('clusterResolution', 0.6)
        
        # --- 检查 igraph 是否可用 ---
        try:
            import igraph
        except ImportError:
            raise ImportError("Please install the igraph package: `conda install -c conda-forge python-igraph` or `pip3 install igraph`.")
        # --- 检查 igraph 是否可用 ---
        
        logger.info(f"Performing clustering using {cluster_method} with resolution={resolution}...")
        if cluster_method == 'leiden':
            sc.tl.leiden(self.adata, resolution=resolution, flavor="igraph", n_iterations=2)
        elif cluster_method == 'louvain':
            sc.tl.louvain(self.adata, resolution=resolution)
        
        logger.info(f"Clustering completed using {cluster_method}, resolution={resolution}")
        self.update_progress(60, 'running')
        return self.adata
    
    def cell_annotation(self, params):
        self.update_progress(65, 'running')
        logger.info("Starting cell annotation...")
        
        try:
            logger.info("Ranking marker genes for clusters...")
            sc.tl.rank_genes_groups(self.adata, 'leiden', method='wilcoxon')
        except Exception as e:
            logger.warning(f"Could not perform rank_genes_groups: {e}")
        
        self.update_progress(75, 'running')
        return self.adata
    
    def differential_expression(self, params):
        self.update_progress(80, 'running')
        logger.info("Starting differential expression analysis...")
        
        try:
            logger.info("Re-ranking genes for DE analysis...")
            sc.tl.rank_genes_groups(self.adata, 'leiden', method='t-test')
        except Exception as e:
            logger.warning(f"Could not perform differential expression: {e}")
        
        self.update_progress(90, 'running')
        return self.adata
    
    def save_results(self):
        logger.info("Saving results...")
        
        # 保存处理后的AnnData
        adata_filename = 'processed_data.h5ad'
        adata_path = os.path.join(self.results_dir, adata_filename)
        logger.info(f"Saving processed AnnData to {adata_path}")
        self.adata.write(adata_path)
        
        # 保存差异表达基因
        de_filename = 'differential_genes.csv'
        de_path = os.path.join(self.results_dir, de_filename)
        logger.info(f"Saving differential genes to {de_path}")
        try:
            de_df = sc.get.rank_genes_groups_df(self.adata, group=None)
            de_df.to_csv(de_path, index=False)
        except Exception as e:
            logger.warning(f"Could not save differential genes: {e}")
            # 创建空的CSV文件
            pd.DataFrame().to_csv(de_path, index=False)
        
        # 保存UMAP可视化
        umap_filename = 'umap_visualization.png'
        umap_path = os.path.join(self.results_dir, umap_filename)
        logger.info(f"Saving UMAP visualization to {umap_path}")
        try:
            plt.figure(figsize=(10, 8))
            sc.pl.umap(self.adata, color='leiden', show=False)
            plt.title("UMAP Clustering (Leiden)")
            plt.tight_layout()
            plt.savefig(umap_path, dpi=300, bbox_inches='tight')
            plt.close()
        except Exception as e:
            logger.error(f"Failed to save UMAP visualization: {e}")
            # 创建一个简单的占位图
            plt.figure(figsize=(10, 8))
            plt.text(0.5, 0.5, 'Visualization Error', horizontalalignment='center', verticalalignment='center')
            plt.title("Visualization Error")
            plt.savefig(umap_path, dpi=300, bbox_inches='tight')
            plt.close()
        
        # 保存聚类结果
        cluster_filename = 'cluster_annotations.csv'
        cluster_path = os.path.join(self.results_dir, cluster_filename)
        logger.info(f"Saving cluster annotations to {cluster_path}")
        try:
            cluster_df = pd.DataFrame({
                'cell_barcode': self.adata.obs.index,
                'cluster': self.adata.obs['leiden'].astype(str),
                'umap_1': self.adata.obsm['X_umap'][:, 0],
                'umap_2': self.adata.obsm['X_umap'][:, 1]
            })
            cluster_df.to_csv(cluster_path, index=False)
        except Exception as e:
            logger.error(f"Failed to save cluster annotations: {e}")
            pd.DataFrame().to_csv(cluster_path, index=False)
        
        # 创建结果对象
        results = []
        
        # 保存AnnData文件
        with open(adata_path, 'rb') as f:
            content = f.read()
            file_path = default_storage.save(
                f'sc_results/{self.task.task_id}/{adata_filename}', 
                ContentFile(content)
            )
            results.append(ScAnalysisResult(
                task=self.task,
                result_type='anndata',
                file=file_path,
                filename=adata_filename
            ))
        
        # 保存差异表达基因文件
        with open(de_path, 'rb') as f:
            content = f.read()
            file_path = default_storage.save(
                f'sc_results/{self.task.task_id}/{de_filename}', 
                ContentFile(content)
            )
            results.append(ScAnalysisResult(
                task=self.task,
                result_type='de-genes',
                file=file_path,
                filename=de_filename
            ))
        
        # 保存UMAP可视化文件
        with open(umap_path, 'rb') as f:
            content = f.read()
            file_path = default_storage.save(
                f'sc_results/{self.task.task_id}/{umap_filename}', 
                ContentFile(content)
            )
            results.append(ScAnalysisResult(
                task=self.task,
                result_type='visualizations',
                file=file_path,
                filename=umap_filename
            ))
        
        # 保存聚类注释文件
        with open(cluster_path, 'rb') as f:
            content = f.read()
            file_path = default_storage.save(
                f'sc_results/{self.task.task_id}/{cluster_filename}', 
                ContentFile(content)
            )
            results.append(ScAnalysisResult(
                task=self.task,
                result_type='cell-annotations',
                file=file_path,
                filename=cluster_filename
            ))
        
        logger.info(f"Creating {len(results)} result entries in database...")
        ScAnalysisResult.objects.bulk_create(results)
        self.update_progress(100, 'completed')
        self.task.completed_at = datetime.now()
        self.task.save()
        
        logger.info(f"Analysis task {self.task.task_id} completed successfully")
        return results


# 视图函数（增强版）
def analysis_page(request):
    """返回单细胞分析页面"""
    logger.info("Serving analysis page.")
    return render(request, 'hormone_app/analysis.html')


@csrf_exempt
def create_task(request):
    """创建新的分析任务"""
    if request.method == 'POST':
        try:
            task_name = request.POST.get('task_name', '未命名分析')
            logger.info(f"Creating new task: {task_name}")
            task = ScAnalysisTask.objects.create(
                user=request.user if request.user.is_authenticated else None,
                task_name=task_name
            )
            logger.info(f"Created new task: {task.task_id}")
            return JsonResponse({
                'id': str(task.task_id), # 添加 'id' 字段
                'task_id': str(task.task_id),
                'status': task.status,
                'created_at': task.created_at.isoformat(),
                'message': '任务创建成功'
            })
        except Exception as e:
            logger.error(f"Error creating task: {e}")
            return JsonResponse({'error': str(e)}, status=500)
    return JsonResponse({'error': '无效请求'}, status=400)


@csrf_exempt
def upload_file(request, task_id):
    """上传分析文件"""
    try:
        task = ScAnalysisTask.objects.get(task_id=task_id)
        logger.info(f"Uploading file for task {task_id}")
    except ScAnalysisTask.DoesNotExist:
        logger.error(f"Task {task_id} does not exist.")
        return JsonResponse({'error': '任务不存在'}, status=404)
    
    if request.method == 'POST' and request.FILES.get('file'):
        try:
            file = request.FILES['file']
            file_type = request.POST.get('file_type', 'sc-data')
            
            # 检查文件大小限制
            max_size = 100 * 1024 * 1024  # 100MB
            if file.size > max_size:
                logger.warning(f"File {file.name} is too large: {file.size} bytes > {max_size} bytes.")
                return JsonResponse({'error': '文件过大，最大支持100MB'}, status=400)
            
            file_path = default_storage.save(
                f'sc_uploads/{task_id}/{file.name}', 
                ContentFile(file.read())
            )
            
            uploaded_file = ScUploadedFile.objects.create(
                task=task,
                file_type=file_type,
                file=file_path,
                filename=file.name,
                file_size=file.size
            )
            
            logger.info(f"File uploaded: {file.name} for task {task_id}")
            return JsonResponse({
                'file_id': str(uploaded_file.file_id),
                'filename': uploaded_file.filename,
                'file_size': uploaded_file.file_size,
                'uploaded_at': uploaded_file.uploaded_at.isoformat(),
                'message': '文件上传成功'
            })
        except Exception as e:
            logger.error(f"Error uploading file: {e}")
            return JsonResponse({'error': str(e)}, status=500)
    
    return JsonResponse({'error': '无效请求'}, status=400)


@csrf_exempt
def run_analysis(request, task_id):
    """运行单细胞分析"""
    try:
        task = ScAnalysisTask.objects.get(task_id=task_id)
        logger.info(f"Starting analysis for task {task_id}")
    except ScAnalysisTask.DoesNotExist:
        logger.error(f"Task {task_id} does not exist.")
        return JsonResponse({'error': '任务不存在'}, status=404)
    
    if request.method == 'POST':
        try:
            task.status = 'running'
            task.progress = 0
            task.save()
            logger.info(f"Set task {task_id} status to running.")
            
            import json
            params = json.loads(request.body) if request.body else {}
            logger.info(f"Received analysis parameters: {params.keys()}")
            
            # 更新任务参数
            for param_type in ['qc', 'normalization', 'clustering', 'annotation', 'de', 'enrichment']:
                if param_type in params:
                    setattr(task, f'{param_type}_params', params[param_type])
            task.save()
            logger.info(f"Updated task {task_id} parameters in database.")
            
            # 查找单细胞数据文件
            sc_data_files = task.uploaded_files.filter(file_type='sc-data')
            if not sc_data_files.exists():
                logger.error("No single-cell data file found for analysis.")
                raise ValueError("未找到单细胞数据文件，请先上传")
            sc_data_file = sc_data_files.first()
            logger.info(f"Found data file: {sc_data_file.filename}")
            
            def _background_sc_analysis():
                try:
                    bg_task = ScAnalysisTask.objects.get(task_id=task_id)
                    analysis_tool = ScAnalysisTool(bg_task)
                    analysis_tool.load_data(sc_data_file.file.path)
                    analysis_tool.quality_control(bg_task.qc_params)
                    analysis_tool.normalization(bg_task.normalization_params)
                    analysis_tool.dimensionality_reduction(bg_task.clustering_params)
                    analysis_tool.cell_annotation(bg_task.annotation_params)
                    analysis_tool.differential_expression(bg_task.de_params)
                    analysis_tool.save_results()
                    logger.info(f"Analysis for task {task_id} finished successfully.")
                except Exception as e:
                    logger.exception(f"Background SC analysis failed for task {task_id}: {e}")
                    ScAnalysisTask.objects.filter(task_id=task_id).update(status='failed', error_message=str(e), progress=0)

            thread = threading.Thread(target=_background_sc_analysis, name=f"SCAnalysis-{task_id}", daemon=True)
            thread.start()
            return JsonResponse({
                'status': 'running',
                'task_id': str(task.task_id),
                'progress': task.progress,
                'message': '分析任务已启动'
            })
            
        except Exception as e:
            logger.error(f"Analysis failed for task {task_id}: {e}")
            task.status = 'failed'
            task.error_message = str(e)
            task.save()
            return JsonResponse({
                'error': str(e),
                'status': 'failed'
            }, status=500)
    
    return JsonResponse({'error': '无效请求'}, status=400)


def get_task_status(request, task_id):
    """获取任务状态"""
    try:
        task = ScAnalysisTask.objects.get(task_id=task_id)
        logger.info(f"Retrieved status for task {task_id}: {task.status}, progress: {task.progress}%")
        return JsonResponse({
            'task_id': str(task.task_id),
            'status': task.status,
            'progress': task.progress,
            'error_message': task.error_message,
            'created_at': task.created_at.isoformat(),
            'completed_at': task.completed_at.isoformat() if task.completed_at else None
        })
    except ScAnalysisTask.DoesNotExist:
        logger.error(f"Task {task_id} does not exist when getting status.")
        return JsonResponse({'error': '任务不存在'}, status=404)


def get_task_results(request, task_id):
    """获取任务结果"""
    try:
        task = ScAnalysisTask.objects.get(task_id=task_id)
        results = task.results.all()
        logger.info(f"Retrieved {len(results)} results for task {task_id}.")
        
        return JsonResponse({
            'task_id': str(task.task_id),
            'status': task.status,
            'total_results': len(results),
            'results': [{
                'result_id': str(res.result_id),
                'result_type': res.result_type,
                'filename': res.filename,
                'download_url': f'/hormone_app/api/results/{res.result_id}/download/',
                'created_at': res.created_at.isoformat()
            } for res in results]
        })
    except ScAnalysisTask.DoesNotExist:
        logger.error(f"Task {task_id} does not exist when getting results.")
        return JsonResponse({'error': '任务不存在'}, status=404)


def download_result(request, result_id):
    """下载分析结果"""
    try:
        result = ScAnalysisResult.objects.get(result_id=result_id)
        logger.info(f"Downloading result file: {result.filename}")
        
        if not os.path.exists(result.file.path):
            logger.error(f"Result file does not exist: {result.file.path}")
            return JsonResponse({'error': '文件不存在'}, status=404)
        
        with open(result.file.path, 'rb') as f:
            response = HttpResponse(f.read(), content_type='application/octet-stream')
            response['Content-Disposition'] = f'attachment; filename="{result.filename}"'
            return response
            
    except ScAnalysisResult.DoesNotExist:
        logger.error(f"Result {result_id} does not exist.")
        return JsonResponse({'error': '结果不存在'}, status=404)


def download_all_results(request, task_id):
    """下载所有分析结果"""
    try:
        task = ScAnalysisTask.objects.get(task_id=task_id)
        results = task.results.all()
        logger.info(f"Downloading all {len(results)} results for task {task_id}.")
        
        if not results:
            logger.warning(f"No results found for task {task_id} to download.")
            return JsonResponse({'error': '没有结果文件'}, status=404)
        
        with tempfile.NamedTemporaryFile(delete=False, suffix='.zip') as tmp_file:
            with zipfile.ZipFile(tmp_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
                for res in results:
                    if os.path.exists(res.file.path):
                        # 使用相对路径，避免在zip中包含完整路径
                        arcname = os.path.basename(res.file.path)
                        zipf.write(res.file.path, arcname)
                        logger.debug(f"Added {arcname} to zip archive.")
            
            with open(tmp_file.name, 'rb') as f:
                response = HttpResponse(f.read(), content_type='application/zip')
                response['Content-Disposition'] = f'attachment; filename="analysis_results_{task_id}.zip"'
                
            os.unlink(tmp_file.name)
            logger.info(f"Downloaded zip archive for task {task_id}.")
            return response
            
    except ScAnalysisTask.DoesNotExist:
        logger.error(f"Task {task_id} does not exist when downloading all results.")
        return JsonResponse({'error': '任务不存在'}, status=404)


# --- 新增的视图函数 ---

def get_visualization_plot(request, task_id, plot_type):
    """生产级可视化接口：返回Plotly JSON；若已有静态PNG则返回image_url。"""
    try:
        task = ScAnalysisTask.objects.get(task_id=task_id)
    except ScAnalysisTask.DoesNotExist:
        return JsonResponse({'error': '任务不存在'}, status=404)
    results_dir = os.path.join(sc_results_dir, str(task_id))
    adata_path = os.path.join(results_dir, 'processed_data.h5ad')
    try:
        # 对大型h5ad读取只发生在用户请求图表时，避免分析主链路阻塞
        if os.path.exists(adata_path):
            adata = sc.read_h5ad(adata_path)
        else:
            res = task.results.filter(result_type='anndata').first()
            if not res or not os.path.exists(res.file.path):
                return JsonResponse({'error': '处理后的AnnData结果不存在'}, status=404)
            adata = sc.read_h5ad(res.file.path)

        if plot_type in ['dim-reduction', 'umap', 'cell-types']:
            if 'X_umap' not in adata.obsm:
                return JsonResponse({'error': 'UMAP坐标不存在'}, status=404)
            clusters = adata.obs['leiden'].astype(str).tolist() if 'leiden' in adata.obs else ['cell'] * adata.n_obs
            data = [{
                'x': adata.obsm['X_umap'][:, 0].astype(float).tolist(),
                'y': adata.obsm['X_umap'][:, 1].astype(float).tolist(),
                'type': 'scattergl', 'mode': 'markers',
                'marker': {'size': 5, 'opacity': 0.75, 'color': clusters},
                'text': adata.obs_names.astype(str).tolist(),
                'name': 'cells'
            }]
            return JsonResponse({'success': True, 'plot_data': {'data': data, 'layout': {'title': 'UMAP Clustering', 'xaxis': {'title': 'UMAP1'}, 'yaxis': {'title': 'UMAP2'}, 'height': 600}}})

        if plot_type == 'qc-metrics':
            obs = adata.obs
            metrics = [c for c in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'] if c in obs]
            data = [{'y': pd.to_numeric(obs[c], errors='coerce').dropna().astype(float).tolist(), 'type': 'box', 'name': c} for c in metrics]
            return JsonResponse({'success': True, 'plot_data': {'data': data, 'layout': {'title': 'Quality Control Metrics', 'height': 500}}})

        if plot_type == 'gene-expression':
            gene = request.GET.get('gene') or (adata.var_names[0] if adata.n_vars else None)
            if not gene or gene not in adata.var_names:
                return JsonResponse({'error': '基因不存在'}, status=404)
            vals = adata[:, gene].X
            vals = vals.toarray().ravel() if hasattr(vals, 'toarray') else np.asarray(vals).ravel()
            data = [{'x': adata.obsm['X_umap'][:,0].astype(float).tolist(), 'y': adata.obsm['X_umap'][:,1].astype(float).tolist(), 'type':'scattergl','mode':'markers','marker':{'size':5,'color':vals.astype(float).tolist(),'colorscale':'Viridis','showscale':True}}]
            return JsonResponse({'success': True, 'plot_data': {'data': data, 'layout': {'title': f'Expression: {gene}', 'height': 600}}})

        if plot_type in ['diff-genes', 'enrichment']:
            de_path = os.path.join(results_dir, 'differential_genes.csv')
            if not os.path.exists(de_path):
                return JsonResponse({'error': '差异基因结果不存在'}, status=404)
            de = pd.read_csv(de_path).head(50)
            if plot_type == 'diff-genes':
                x = de['logfoldchanges'].astype(float).tolist() if 'logfoldchanges' in de else []
                pcol = 'pvals_adj' if 'pvals_adj' in de else ('pvals' if 'pvals' in de else None)
                y = (-np.log10(pd.to_numeric(de[pcol], errors='coerce').fillna(1).clip(lower=1e-300))).astype(float).tolist() if pcol else []
                return JsonResponse({'success': True, 'plot_data': {'data': [{'x': x, 'y': y, 'type': 'scatter', 'mode': 'markers', 'text': de.get('names', pd.Series(range(len(de)))).astype(str).tolist()}], 'layout': {'title': 'Differential Expression Volcano', 'xaxis': {'title': 'logFC'}, 'yaxis': {'title': '-log10(FDR/p)'}}}})
            names = de.get('names', pd.Series([])).astype(str).head(20).tolist()
            scores = pd.to_numeric(de.get('scores', pd.Series([0]*len(names))), errors='coerce').fillna(0).head(20).astype(float).tolist()
            return JsonResponse({'success': True, 'plot_data': {'data': [{'x': scores, 'y': names, 'type': 'bar', 'orientation': 'h'}], 'layout': {'title': 'Top Marker Genes', 'height': 600}}})

        # 静态图兜底：只返回真实存在的文件，拒绝占位图
        filename = {'heatmap':'gene_heatmap.png', 'cell-composition':'cell_composition.png', 'trajectory':'trajectory.png'}.get(plot_type)
        if filename and os.path.exists(os.path.join(results_dir, filename)):
            return JsonResponse({'success': True, 'plot_data': {'image_url': f'/media/sc_results/{task_id}/{filename}'}})
        return JsonResponse({'error': f'图表暂不可用或未生成: {plot_type}'}, status=404)
    except Exception as e:
        logger.exception(f"Visualization error for task {task_id}, plot {plot_type}: {e}")
        return JsonResponse({'error': str(e)}, status=500)


def download_result_by_type(request, task_id):
    """
    下载指定任务的特定类型结果。
    从前端的URL: /hormone_app/api/results/<task_id>/download/?type=<result_type>
    获取result_type参数，在数据库中查找对应的结果记录，然后调用现有的下载逻辑。
    """
    try:
        task = ScAnalysisTask.objects.get(task_id=task_id)
        logger.info(f"Downloading result by type for task {task_id}, type: {request.GET.get('type')}.")
    except ScAnalysisTask.DoesNotExist:
        logger.error(f"Task {task_id} does not exist when downloading result by type.")
        return JsonResponse({'error': '任务不存在'}, status=404)

    result_type = request.GET.get('type')
    if not result_type:
        logger.warning("Missing 'type' parameter when downloading result by type.")
        return JsonResponse({'error': '缺少必需的 "type" 参数'}, status=400)

    # 检查result_type是否有效
    valid_types = [choice[0] for choice in ScAnalysisResult.RESULT_TYPES]
    if result_type not in valid_types:
        logger.warning(f"Invalid result type requested: {result_type}")
        return JsonResponse({'error': f'无效的结果类型: {result_type}'}, status=400)

    # 在指定任务的结果中查找第一个匹配的类型
    try:
        result = task.results.get(result_type=result_type)
        logger.info(f"Found result {result.result_id} of type {result_type} for task {task_id}.")
        # 调用现有的 download_result 函数处理下载
        return download_result(request, result.result_id)
    except ScAnalysisResult.DoesNotExist:
        logger.error(f"No result of type '{result_type}' found for task {task_id}.")
        return JsonResponse({'error': f'任务 {task_id} 中找不到类型为 "{result_type}" 的结果'}, status=404)
    except ScAnalysisResult.MultipleObjectsReturned:
        # 如果有多个同类型的文件，可以选择下载第一个或返回列表让用户选择
        # 此处我们下载第一个
        logger.warning(f"Multiple results of type '{result_type}' found for task {task_id}. Downloading first one.")
        result = task.results.filter(result_type=result_type).first()
        return download_result(request, result.result_id)


# 基因表达分析相关（示例结构）
class GeneExpressionAnalysis:
    def gene_analysis_page(self, request):
        """基因表达分析页面"""
        logger.info("Serving gene expression analysis page.")
        return render(request, 'hormone_app/gene_expression_analysis.html')

    def create_gene_task(self, request):
        """创建基因表达分析任务"""
        if request.method == 'POST':
            try:
                task_name = request.POST.get('task_name', '未命名基因表达分析')
                logger.info(f"Creating new gene expression task: {task_name}")
                task = ScAnalysisTask.objects.create(
                    user=request.user if request.user.is_authenticated else None,
                    task_name=task_name
                )
                return JsonResponse({
                    'task_id': str(task.task_id),
                    'status': task.status,
                    'created_at': task.created_at.isoformat(),
                    'message': '基因表达分析任务创建成功'
                })
            except Exception as e:
                logger.error(f"Error creating gene expression task: {e}")
                return JsonResponse({'error': str(e)}, status=500)
        return JsonResponse({'error': '无效请求'}, status=400)

    def upload_gene_file(self, request, task_id):
        """上传基因表达分析文件"""
        try:
            task = ScAnalysisTask.objects.get(task_id=task_id)
            logger.info(f"Uploading gene expression file for task {task_id}")
        except ScAnalysisTask.DoesNotExist:
            logger.error(f"Task {task_id} does not exist for gene file upload.")
            return JsonResponse({'error': '任务不存在'}, status=404)
        
        if request.method == 'POST' and request.FILES.get('file'):
            try:
                file = request.FILES['file']
                file_type = request.POST.get('file_type', 'sc-data')
                
                file_path = default_storage.save(
                    f'sc_uploads/{task_id}/{file.name}', 
                    ContentFile(file.read())
                )
                
                uploaded_file = ScUploadedFile.objects.create(
                    task=task,
                    file_type=file_type,
                    file=file_path,
                    filename=file.name,
                    file_size=file.size
                )
                
                logger.info(f"Gene file uploaded: {file.name} for task {task_id}")
                return JsonResponse({
                    'file_id': str(uploaded_file.file_id),
                    'filename': uploaded_file.filename,
                    'uploaded_at': uploaded_file.uploaded_at.isoformat()
                })
            except Exception as e:
                logger.error(f"Error uploading gene expression file: {e}")
                return JsonResponse({'error': str(e)}, status=500)
        
        return JsonResponse({'error': '无效请求'}, status=400)

    def run_gene_analysis(self, request, task_id):
        """运行基因表达分析"""
        # 这里可以实现特定的基因表达分析逻辑
        logger.info(f"Running gene expression analysis for task {task_id} (stub implementation).")
        return JsonResponse({'message': '基因表达分析功能待实现'})

    def get_gene_task_status(self, request, task_id):
        """获取基因表达分析任务状态"""
        logger.info(f"Getting gene expression task status for {task_id}.")
        return get_task_status(request, task_id)

    def get_gene_task_results(self, request, task_id):
        """获取基因表达分析结果"""
        logger.info(f"Getting gene expression task results for {task_id}.")
        return get_task_results(request, task_id)

    def download_gene_result(self, request, result_id):
        """下载基因表达分析结果"""
        logger.info(f"Downloading gene expression result {result_id}.")
        return download_result(request, result_id)

    def download_all_gene_results(self, request, task_id):
        """下载所有基因表达分析结果"""
        logger.info(f"Downloading all gene expression results for task {task_id}.")
        return download_all_results(request, task_id)

# 实例化基因表达分析类
gene_expression_analysis = GeneExpressionAnalysis()