import os
import uuid
import tempfile
import zipfile
from datetime import datetime
from django.http import JsonResponse, HttpResponse
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from django.conf import settings
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from django.contrib.auth.models import User
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

# 从models.py导入所需模型（解决冲突的关键）
from ..models import GeneAnalysisTask, GeneUploadedFile, GeneAnalysisResult

# 确保媒体文件目录存在
os.makedirs(os.path.join(settings.MEDIA_ROOT, 'gene_uploads'), exist_ok=True)
os.makedirs(os.path.join(settings.MEDIA_ROOT, 'gene_results'), exist_ok=True)

# 分析工具类
class GeneAnalysisTool:
    def __init__(self, task):
        self.task = task
        self.expression_data = None
        self.sample_info = None
        self.gene_annotations = None
        self.results_dir = os.path.join(settings.MEDIA_ROOT, 'gene_results', str(task.task_id))
        os.makedirs(self.results_dir, exist_ok=True)
        
        # 设置matplotlib中文字体支持
        plt.rcParams["font.family"] = ["SimHei", "WenQuanYi Micro Hei", "Heiti TC"]
        plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题
    
    def update_progress(self, progress, status_text=None):
        """更新任务进度"""
        self.task.progress = progress
        if status_text:
            self.task.status = status_text
        self.task.save()
    
    def load_data(self):
        """加载分析所需数据"""
        self.update_progress(10, 'running')
        
        # 获取表达数据文件
        expr_file = self.task.uploaded_files.filter(file_type='expression-data').first()
        if not expr_file:
            raise ValueError("未找到表达数据文件")
        
        # 读取表达数据
        file_ext = os.path.splitext(expr_file.filename)[1].lower()
        if file_ext == '.csv':
            self.expression_data = pd.read_csv(expr_file.file.path, index_col=0)
        elif file_ext == '.tsv' or file_ext == '.txt':
            self.expression_data = pd.read_csv(expr_file.file.path, sep='\t', index_col=0)
        elif file_ext == '.xlsx':
            self.expression_data = pd.read_excel(expr_file.file.path, index_col=0)
        else:
            raise ValueError(f"不支持的表达数据格式: {file_ext}")
        
        # 获取样本信息（如果有）
        sample_file = self.task.uploaded_files.filter(file_type='sample-info').first()
        if sample_file:
            file_ext = os.path.splitext(sample_file.filename)[1].lower()
            if file_ext == '.csv':
                self.sample_info = pd.read_csv(sample_file.file.path, index_col=0)
            elif file_ext == '.tsv' or file_ext == '.txt':
                self.sample_info = pd.read_csv(sample_file.file.path, sep='\t', index_col=0)
            elif file_ext == '.xlsx':
                self.sample_info = pd.read_excel(sample_file.file.path, index_col=0)
        
        self.update_progress(20, 'running')
        return True
    
    def normalize_data(self):
        """数据标准化处理"""
        self.update_progress(25, 'running')
        
        method = self.task.normalization_params.get('method', 'tpm')
        
        if method == 'tpm':
            # TPM标准化
            row_sums = self.expression_data.sum(axis=1)
            self.expression_data = self.expression_data.div(row_sums, axis=0) * 1e6
            self.expression_data = np.log2(self.expression_data + 1)  # 加1避免log(0)
            
        elif method == 'zscore':
            # Z-score标准化
            self.expression_data = (self.expression_data - self.expression_data.mean()) / self.expression_data.std()
            
        elif method == 'quantile':
            # 分位数标准化
            rank_mean = self.expression_data.stack().groupby(self.expression_data.rank(method='first').stack().astype(int)).mean()
            self.expression_data = self.expression_data.rank(method='min').stack().astype(int).map(rank_mean).unstack()
        
        self.update_progress(35, 'running')
        return True
    
    def differential_expression(self):
        """差异表达分析"""
        self.update_progress(40, 'running')
        
        # 获取分组信息
        group_col = self.task.de_params.get('group_column', 'group')
        group1 = self.task.de_params.get('group1', 'control')
        group2 = self.task.de_params.get('group2', 'treatment')
        
        if not self.sample_info or group_col not in self.sample_info.columns:
            raise ValueError(f"样本信息中未找到分组列: {group_col}")
        
        # 获取两组样本
        samples1 = self.sample_info[self.sample_info[group_col] == group1].index.tolist()
        samples2 = self.sample_info[self.sample_info[group_col] == group2].index.tolist()
        
        # 确保样本存在于表达数据中
        samples1 = [s for s in samples1 if s in self.expression_data.columns]
        samples2 = [s for s in samples2 if s in self.expression_data.columns]
        
        if len(samples1) == 0 or len(samples2) == 0:
            raise ValueError("未找到有效的样本数据用于差异分析")
        
        # 执行差异分析
        de_results = []
        for gene, row in self.expression_data.iterrows():
            # 获取两组的表达值
            expr1 = row[samples1].values
            expr2 = row[samples2].values
            
            # 进行t检验
            stat, p_value = stats.ttest_ind(expr1, expr2, equal_var=False)
            
            # 计算平均值和倍数变化
            mean1 = np.mean(expr1)
            mean2 = np.mean(expr2)
            log2fc = mean2 - mean1  # 已经过log转换
            
            de_results.append({
                'gene': gene,
                'mean_group1': mean1,
                'mean_group2': mean2,
                'log2fc': log2fc,
                'p_value': p_value,
                'significant': 'yes' if p_value < 0.05 and abs(log2fc) > 1 else 'no'
            })
        
        # 转换为DataFrame并校正p值
        self.de_df = pd.DataFrame(de_results)
        self.de_df['adj_p_value'] = stats.false_discovery_control(self.de_df['p_value'])
        
        self.update_progress(60, 'running')
        return True
    
    def generate_visualizations(self):
        """生成可视化结果"""
        self.update_progress(65, 'running')
        
        # 1. 火山图
        plt.figure(figsize=(10, 8))
        sig = self.de_df[(self.de_df['adj_p_value'] < 0.05) & (abs(self.de_df['log2fc']) > 1)]
        non_sig = self.de_df[(self.de_df['adj_p_value'] >= 0.05) | (abs(self.de_df['log2fc']) <= 1)]
        
        plt.scatter(non_sig['log2fc'], -np.log10(non_sig['adj_p_value']), 
                   color='gray', alpha=0.5, label='不显著')
        plt.scatter(sig['log2fc'], -np.log10(sig['adj_p_value']), 
                   color='red', alpha=0.7, label='显著')
        
        plt.axhline(y=-np.log10(0.05), color='black', linestyle='--')
        plt.axvline(x=1, color='black', linestyle='--')
        plt.axvline(x=-1, color='black', linestyle='--')
        
        plt.xlabel('log2(倍数变化)')
        plt.ylabel('-log10(校正后p值)')
        plt.title('差异表达基因火山图')
        plt.legend()
        plt.tight_layout()
        
        volcano_path = os.path.join(self.results_dir, 'volcano_plot.png')
        plt.savefig(volcano_path, dpi=300)
        plt.close()
        
        # 2. 热图（使用前50个显著差异基因）
        top_genes = self.de_df.sort_values('adj_p_value').head(50)['gene'].tolist()
        if len(top_genes) >= 5:  # 确保有足够的基因绘制热图
            plt.figure(figsize=(12, 10))
            heatmap_data = self.expression_data.loc[top_genes]
            
            # 标准化行（基因）
            heatmap_data = (heatmap_data - heatmap_data.mean(axis=1)[:, np.newaxis]) / heatmap_data.std(axis=1)[:, np.newaxis]
            
            sns.heatmap(heatmap_data, cmap='coolwarm', center=0)
            plt.title('Top 50 差异表达基因热图')
            plt.tight_layout()
            
            heatmap_path = os.path.join(self.results_dir, 'heatmap.png')
            plt.savefig(heatmap_path, dpi=300)
            plt.close()
        else:
            heatmap_path = None
        
        self.update_progress(80, 'running')
        return {
            'volcano_plot': volcano_path,
            'heatmap': heatmap_path
        }
    
    def save_results(self, visualizations):
        """保存分析结果"""
        results = []
        
        # 保存标准化数据
        norm_path = os.path.join(self.results_dir, 'normalized_expression.csv')
        self.expression_data.to_csv(norm_path)
        
        with open(norm_path, 'rb') as f:
            content = f.read()
            file_name = os.path.basename(norm_path)
            file_path = default_storage.save(f'gene_results/{self.task.task_id}/{file_name}', ContentFile(content))
            results.append(GeneAnalysisResult(
                task=self.task,
                result_type='normalized-data',
                file=file_path,
                filename=file_name
            ))
        
        # 保存差异分析结果
        de_path = os.path.join(self.results_dir, 'differential_expression_results.csv')
        self.de_df.to_csv(de_path, index=False)
        
        with open(de_path, 'rb') as f:
            content = f.read()
            file_name = os.path.basename(de_path)
            file_path = default_storage.save(f'gene_results/{self.task.task_id}/{file_name}', ContentFile(content))
            results.append(GeneAnalysisResult(
                task=self.task,
                result_type='de-results',
                file=file_path,
                filename=file_name
            ))
        
        # 保存火山图
        if 'volcano_plot' in visualizations and visualizations['volcano_plot']:
            with open(visualizations['volcano_plot'], 'rb') as f:
                content = f.read()
                file_name = os.path.basename(visualizations['volcano_plot'])
                file_path = default_storage.save(f'gene_results/{self.task.task_id}/{file_name}', ContentFile(content))
                results.append(GeneAnalysisResult(
                    task=self.task,
                    result_type='volcano-plot',
                    file=file_path,
                    filename=file_name
                ))
        
        # 保存热图
        if 'heatmap' in visualizations and visualizations['heatmap']:
            with open(visualizations['heatmap'], 'rb') as f:
                content = f.read()
                file_name = os.path.basename(visualizations['heatmap'])
                file_path = default_storage.save(f'gene_results/{self.task.task_id}/{file_name}', ContentFile(content))
                results.append(GeneAnalysisResult(
                    task=self.task,
                    result_type='heatmap',
                    file=file_path,
                    filename=file_name
                ))
        
        # 批量保存结果记录
        GeneAnalysisResult.objects.bulk_create(results)
        
        self.update_progress(100, 'completed')
        self.task.completed_at = datetime.now()
        self.task.save()
        
        return results

# 视图函数
def gene_analysis_page(request):
    """基因分析页面视图"""
    return render(request, 'hormone_app/gene_analysis.html')

@csrf_exempt
def create_gene_task(request):
    """创建新的基因分析任务"""
    if request.method == 'POST':
        task_name = request.POST.get('task_name', '未命名基因分析')
        task = GeneAnalysisTask.objects.create(
            user=request.user if request.user.is_authenticated else None,
            task_name=task_name
        )
        return JsonResponse({
            'task_id': str(task.task_id),
            'status': task.status,
            'created_at': task.created_at.isoformat()
        })
    return JsonResponse({'error': '无效请求'}, status=400)

@csrf_exempt
def upload_gene_file(request, task_id):
    """上传文件到指定基因分析任务"""
    try:
        task = GeneAnalysisTask.objects.get(task_id=task_id)
    except GeneAnalysisTask.DoesNotExist:
        return JsonResponse({'error': '任务不存在'}, status=404)
    
    if request.method == 'POST' and request.FILES.get('file'):
        file = request.FILES['file']
        file_type = request.POST.get('file_type', 'expression-data')
        
        # 保存文件
        file_path = default_storage.save(f'gene_uploads/{task_id}/{file.name}', ContentFile(file.read()))
        
        # 创建文件记录
        uploaded_file = GeneUploadedFile.objects.create(
            task=task,
            file_type=file_type,
            file=file_path,
            filename=file.name,
            file_size=file.size
        )
        
        return JsonResponse({
            'file_id': str(uploaded_file.file_id),
            'filename': uploaded_file.filename,
            'uploaded_at': uploaded_file.uploaded_at.isoformat()
        })
    
    return JsonResponse({'error': '无效请求'}, status=400)

@csrf_exempt
def run_gene_analysis(request, task_id):
    """运行基因分析任务"""
    try:
        task = GeneAnalysisTask.objects.get(task_id=task_id)
    except GeneAnalysisTask.DoesNotExist:
        return JsonResponse({'error': '任务不存在'}, status=404)
    
    if request.method == 'POST':
        # 更新任务状态
        task.status = 'running'
        task.progress = 0
        task.save()
        
        try:
            # 获取分析参数
            import json
            params = json.loads(request.body)
            
            # 保存参数
            if 'data_params' in params:
                task.data_params = params['data_params']
            if 'normalization_params' in params:
                task.normalization_params = params['normalization_params']
            if 'de_params' in params:
                task.de_params = params['de_params']
            if 'visualization_params' in params:
                task.visualization_params = params['visualization_params']
            task.save()
            
            # 执行分析
            analysis_tool = GeneAnalysisTool(task)
            analysis_tool.load_data()
            analysis_tool.normalize_data()
            analysis_tool.differential_expression()
            visualizations = analysis_tool.generate_visualizations()
            analysis_tool.save_results(visualizations)
            
            return JsonResponse({
                'status': 'success',
                'task_id': str(task.task_id),
                'progress': task.progress
            })
            
        except Exception as e:
            # 记录错误信息
            task.status = 'failed'
            task.error_message = str(e)
            task.save()
            return JsonResponse({
                'error': str(e),
                'status': 'failed'
            }, status=500)
    
    return JsonResponse({'error': '无效请求'}, status=400)

def get_gene_task_status(request, task_id):
    """获取基因分析任务状态"""
    try:
        task = GeneAnalysisTask.objects.get(task_id=task_id)
        return JsonResponse({
            'task_id': str(task.task_id),
            'status': task.status,
            'progress': task.progress,
            'error_message': task.error_message,
            'completed_at': task.completed_at.isoformat() if task.completed_at else None
        })
    except GeneAnalysisTask.DoesNotExist:
        return JsonResponse({'error': '任务不存在'}, status=404)

def get_gene_task_results(request, task_id):
    """获取基因分析任务结果列表"""
    try:
        task = GeneAnalysisTask.objects.get(task_id=task_id)
        results = task.results.all()
        
        return JsonResponse({
            'task_id': str(task.task_id),
            'status': task.status,
            'results': [{
                'result_id': str(res.result_id),
                'result_type': res.result_type,
                'filename': res.filename,
                'download_url': f'/gene-analysis/download/{res.result_id}/',
                'created_at': res.created_at.isoformat()
            } for res in results]
        })
    except GeneAnalysisTask.DoesNotExist:
        return JsonResponse({'error': '任务不存在'}, status=404)

def download_gene_result(request, result_id):
    """下载单个基因分析结果文件"""
    try:
        result = GeneAnalysisResult.objects.get(result_id=result_id)
        
        # 检查文件是否存在
        if not os.path.exists(result.file.path):
            return JsonResponse({'error': '文件不存在'}, status=404)
        
        # 读取文件内容
        with open(result.file.path, 'rb') as f:
            response = HttpResponse(f.read(), content_type='application/octet-stream')
            response['Content-Disposition'] = f'attachment; filename="{result.filename}"'
            return response
            
    except GeneAnalysisResult.DoesNotExist:
        return JsonResponse({'error': '结果不存在'}, status=404)

def download_all_gene_results(request, task_id):
    """下载所有基因分析结果文件（打包为ZIP）"""
    try:
        task = GeneAnalysisTask.objects.get(task_id=task_id)
        results = task.results.all()
        
        if not results:
            return JsonResponse({'error': '没有结果文件'}, status=404)
        
        # 创建临时ZIP文件
        with tempfile.NamedTemporaryFile(delete=False, suffix='.zip') as tmp_file:
            with zipfile.ZipFile(tmp_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
                for res in results:
                    if os.path.exists(res.file.path):
                        zipf.write(res.file.path, res.filename)
            
            # 读取ZIP文件内容
            with open(tmp_file.name, 'rb') as f:
                response = HttpResponse(f.read(), content_type='application/zip')
                response['Content-Disposition'] = f'attachment; filename="gene_analysis_results_{task_id}.zip"'
                
            # 删除临时文件
            os.unlink(tmp_file.name)
            
            return response
            
    except GeneAnalysisTask.DoesNotExist:
        return JsonResponse({'error': '任务不存在'}, status=404)

# URL配置（需要添加到项目的urls.py中）
"""
在你的项目urls.py中添加:

from django.urls import path
from .analysishou import gene_expression_analysis  # 注意修改导入路径

urlpatterns = [
    # ...其他URL
    path('gene-analysis/', gene_expression_analysis.gene_analysis_page, name='gene_analysis_page'),
    path('api/gene-tasks/', gene_expression_analysis.create_gene_task, name='create_gene_task'),
    path('api/gene-tasks/<uuid:task_id>/upload/', gene_expression_analysis.upload_gene_file, name='upload_gene_file'),
    path('api/gene-tasks/<uuid:task_id>/run/', gene_expression_analysis.run_gene_analysis, name='run_gene_analysis'),
    path('api/gene-tasks/<uuid:task_id>/status/', gene_expression_analysis.get_gene_task_status, name='get_gene_task_status'),
    path('api/gene-tasks/<uuid:task_id>/results/', gene_expression_analysis.get_gene_task_results, name='get_gene_task_results'),
    path('api/gene-results/<uuid:result_id>/download/', gene_expression_analysis.download_gene_result, name='download_gene_result'),
    path('api/gene-tasks/<uuid:task_id>/download-all/', gene_expression_analysis.download_all_gene_results, name='download_all_gene_results'),
]
"""
    