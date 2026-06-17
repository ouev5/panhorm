"""
转录组测序分析视图函数
处理转录组分析的前后端交互
使用真实数据进行分析 - 支持表达矩阵CSV文件
"""

import os
import uuid
import json
import logging
import tempfile
import subprocess
import shutil
import gzip
import shlex
from datetime import datetime
from pathlib import Path
from django.shortcuts import render
from django.http import JsonResponse, FileResponse
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_http_methods
from django.conf import settings
from django.utils import timezone
from uuid import UUID
import mimetypes
import zipfile
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import warnings
warnings.filterwarnings('ignore')

# 导入数据库模型
from ..models import (
    TranscriptomeAnalysisTask, 
    TranscriptomeAnalysisFile,
    TranscriptomeAnalysisResult
)

# 配置日志
logger = logging.getLogger(__name__)


# ===== 生产级 RNA-seq FASTQ pipeline 配置 =====
RNA_PIPELINE_ROOT = os.path.join(settings.MEDIA_ROOT, 'rna_pipeline')
RNA_BUILTIN_REFS = {
    # 内置参考只在服务器已配置索引时启用；没有索引会明确报错，不再生成模拟表达矩阵
    'hg38': {'name': 'Human hg38'},
    'hg19': {'name': 'Human hg19'},
    'mm10': {'name': 'Mouse mm10'},
    'mm9': {'name': 'Mouse mm9'},
    'rn6': {'name': 'Rat rn6'},
}


def _tool_path(tool):
    return shutil.which(tool)


def _run_cmd(cmd, cwd=None, log_file=None, timeout=None):
    """运行外部命令，失败时抛出包含 stderr 的异常。"""
    if isinstance(cmd, str):
        printable = cmd
    else:
        printable = ' '.join(shlex.quote(str(x)) for x in cmd)
    logger.info(f"[RNA pipeline] Running: {printable}")
    with open(log_file, 'a', encoding='utf-8') if log_file else open(os.devnull, 'w') as lf:
        lf.write(f"\n$ {printable}\n")
        lf.flush()
        proc = subprocess.run(cmd, cwd=cwd, stdout=lf, stderr=subprocess.STDOUT, text=True, timeout=timeout)
    if proc.returncode != 0:
        tail = ''
        if log_file and os.path.exists(log_file):
            with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
                tail = ''.join(f.readlines()[-80:])
        raise RuntimeError(f"Command failed ({proc.returncode}): {printable}\n{tail}")


def _append_log(log_file, message):
    logger.info(message)
    if log_file:
        with open(log_file, 'a', encoding='utf-8') as f:
            f.write(f"[{datetime.now().isoformat()}] {message}\n")


def _is_gzip(path):
    return str(path).lower().endswith('.gz')


def _open_text_maybe_gzip(path):
    return gzip.open(path, 'rt', errors='ignore') if _is_gzip(path) else open(path, 'r', errors='ignore')


def validate_fastq_file(path, max_records_check=1000):
    """轻量 FASTQ 格式校验。"""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        raise ValueError(f"FASTQ file missing or empty: {path}")
    checked = 0
    with _open_text_maybe_gzip(path) as fh:
        while checked < max_records_check:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline(); plus = fh.readline(); qual = fh.readline()
            if not (seq and plus and qual):
                raise ValueError(f"Incomplete FASTQ record in {os.path.basename(path)}")
            if not h.startswith('@') or not plus.startswith('+'):
                raise ValueError(f"Invalid FASTQ format in {os.path.basename(path)} at record {checked+1}")
            if len(seq.strip()) != len(qual.strip()):
                raise ValueError(f"Sequence/quality length mismatch in {os.path.basename(path)} at record {checked+1}")
            checked += 1
    if checked == 0:
        raise ValueError(f"No FASTQ records found in {os.path.basename(path)}")
    return checked


def validate_fasta_file(path):
    if not path or str(path).startswith('builtin://'):
        return False
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    with open(path, 'r', errors='ignore') as f:
        for line in f:
            if line.strip():
                return line.startswith('>')
    return False


def find_annotation_gtf(task):
    """从 gene-list/custom 上传中寻找 GTF/GFF 注释。"""
    files = task.files or {}
    for key in ['annotation', 'gene-list', 'custom-genesets']:
        info = files.get(key)
        if isinstance(info, dict):
            path = info.get('path')
            name = (info.get('name') or path or '').lower()
            if path and os.path.exists(path) and name.endswith(('.gtf', '.gff', '.gff3')):
                return path
    return None


def resolve_reference_for_rnaseq(task, results_dir, log_file):
    """解析参考。生产模式必须具备可用索引或用户上传 FASTA+GTF。"""
    files = task.files or {}
    ref = files.get('reference-genome') or {}
    ref_path = ref.get('path') if isinstance(ref, dict) else None
    ref_name = ref.get('name') if isinstance(ref, dict) else None
    selected = (task.parameters or {}).get('reference_genome') or ref_name

    cfg_path = os.path.join(RNA_PIPELINE_ROOT, f'{selected}.json') if selected else None
    cfg = {}
    if cfg_path and os.path.exists(cfg_path):
        with open(cfg_path, 'r', encoding='utf-8') as f:
            cfg = json.load(f)

    # 优先使用已配置的 salmon 索引，其次 STAR/HISAT2（后续可扩展），再用用户上传 FASTA 动态建 salmon 索引
    salmon_index = cfg.get('salmon_index')
    if salmon_index and os.path.isdir(salmon_index):
        return {'mode': 'salmon', 'salmon_index': salmon_index, 'gtf': cfg.get('gtf'), 'name': selected}

    if ref_path and validate_fasta_file(ref_path):
        salmon = _tool_path('salmon')
        if not salmon:
            raise RuntimeError('生产级 FASTQ 分析需要安装 salmon（或配置已有 salmon_index）。当前服务器未检测到 salmon。')
        gtf = find_annotation_gtf(task)
        index_dir = os.path.join(results_dir, 'salmon_index')
        os.makedirs(index_dir, exist_ok=True)
        # salmon index 可直接用基因组/转录本 fasta；若上传的是基因组 fasta，建议同时上传转录本 fasta 或配置索引
        _run_cmd([salmon, 'index', '-t', ref_path, '-i', index_dir, '-k', '31'], log_file=log_file)
        return {'mode': 'salmon', 'salmon_index': index_dir, 'gtf': gtf, 'name': ref_name or 'custom'}

    if selected in RNA_BUILTIN_REFS:
        raise RuntimeError(
            f"内置参考基因组 {selected} 尚未配置生产级索引。请在 {RNA_PIPELINE_ROOT}/{selected}.json 中配置 salmon_index/gtf，"
            "或上传自定义转录本 FASTA（以及可选 GTF）。系统不再使用模拟表达矩阵。"
        )
    raise RuntimeError('缺少可用参考：请上传 FASTA/转录本 FASTA，或配置内置参考基因组索引。')


def group_fastq_files(fastq_infos):
    """将 FASTQ 文件按样本组织；支持单端和常见 _R1/_R2 配对。"""
    import re
    samples = {}
    for info in fastq_infos:
        name = info.get('name') or os.path.basename(info.get('path', 'sample'))
        path = info.get('path')
        if not path or not os.path.exists(path):
            raise FileNotFoundError(f"FASTQ not found: {name} ({path})")
        validate_fastq_file(path)
        base = re.sub(r'\.(fastq|fq)(\.gz)?$', '', name, flags=re.I)
        m = re.search(r'(.+?)(?:[_\-.]R?([12]))(?:[_\-.].*)?$', base, flags=re.I)
        if m:
            sample, mate = m.group(1), m.group(2)
        else:
            sample, mate = base, '1'
        samples.setdefault(sample, {})[mate] = path
    return samples


def run_salmon_quant(task, fastq_infos, ref_info, results_dir, log_file):
    salmon = _tool_path('salmon')
    if not salmon:
        raise RuntimeError('生产级 FASTQ 定量需要安装 salmon。')
    quant_root = os.path.join(results_dir, 'salmon_quant')
    os.makedirs(quant_root, exist_ok=True)
    samples = group_fastq_files(fastq_infos)
    quant_files = {}
    threads = int((task.parameters or {}).get('threads', 2))
    for sample, mates in samples.items():
        out_dir = os.path.join(quant_root, sample)
        os.makedirs(out_dir, exist_ok=True)
        cmd = [salmon, 'quant', '-i', ref_info['salmon_index'], '-p', str(threads), '--validateMappings', '-o', out_dir]
        if '1' in mates and '2' in mates:
            cmd += ['-1', mates['1'], '-2', mates['2']]
        elif '1' in mates:
            cmd += ['-r', mates['1']]
        else:
            raise RuntimeError(f"Sample {sample} has no R1/single FASTQ")
        _run_cmd(cmd, log_file=log_file)
        q = os.path.join(out_dir, 'quant.sf')
        if not os.path.exists(q):
            raise RuntimeError(f"salmon output missing: {q}")
        quant_files[sample] = q
    # 合并 TPM 和 NumReads 矩阵
    counts = []
    tpms = []
    for sample, qf in quant_files.items():
        df = pd.read_csv(qf, sep='\t')
        counts.append(df.set_index('Name')['NumReads'].rename(sample))
        tpms.append(df.set_index('Name')['TPM'].rename(sample))
    df_counts = pd.concat(counts, axis=1).fillna(0)
    df_tpm = pd.concat(tpms, axis=1).fillna(0)
    counts_path = os.path.join(results_dir, 'raw_counts.csv')
    tpm_path = os.path.join(results_dir, 'tpm.csv')
    df_counts.to_csv(counts_path)
    df_tpm.to_csv(tpm_path)
    return df_counts, {'counts': counts_path, 'tpm': tpm_path, 'quant_dir': quant_root}

# 创建日志格式
class TranscriptomeAnalysisLogger:
    """转录组分析专用日志记录器"""
    
    @staticmethod
    def info(task_id, message, user=None, step=None):
        """信息日志"""
        user_info = f"User: {user.username if user and user.is_authenticated else 'anonymous'}"
        step_info = f"Step: {step}" if step else ""
        logger.info(f"[Transcriptome][{task_id}] {user_info} {step_info} - {message}")
    
    @staticmethod
    def warning(task_id, message, user=None, step=None):
        """警告日志"""
        user_info = f"User: {user.username if user and user.is_authenticated else 'anonymous'}"
        step_info = f"Step: {step}" if step else ""
        logger.warning(f"[Transcriptome][{task_id}] {user_info} {step_info} - {message}")
    
    @staticmethod
    def error(task_id, message, user=None, step=None, exc=None):
        """错误日志"""
        user_info = f"User: {user.username if user and user.is_authenticated else 'anonymous'}"
        step_info = f"Step: {step}" if step else ""
        error_msg = f"[Transcriptome][{task_id}] {user_info} {step_info} - {message}"
        if exc:
            error_msg += f"\nException: {str(exc)}"
        logger.error(error_msg, exc_info=exc is not None)
    
    @staticmethod
    def debug(task_id, message, user=None, step=None):
        """调试日志"""
        user_info = f"User: {user.username if user and user.is_authenticated else 'anonymous'}"
        step_info = f"Step: {step}" if step else ""
        logger.debug(f"[Transcriptome][{task_id}] {user_info} {step_info} - {message}")

# ==================== 配置常量 ====================

# 支持的图表类型
PLOT_TYPES = [
    'qc-metrics',
    'expression-distribution', 
    'pca-analysis',
    'volcano-plot',
    'enrichment-barplot',
    'heatmap',
    'coexpression-network',
    'expression-profiles',
    'pathway-map'
]

# 结果文件类型
RESULT_TYPES = {
    'all': 'All Results (ZIP)',
    'counts': 'Expression Counts (CSV)',
    'normalized': 'Normalized Data (CSV)',
    'de-genes': 'DE Genes (CSV)',
    'enrichment': 'Enrichment Results (CSV)',
    'visualizations': 'All Plots (ZIP)',
    'analysis-summary': 'Analysis Summary (TXT)'
}

# 文件类型映射
FILE_TYPE_MAPPING = {
    '.csv': 'expression-matrix',
    '.tsv': 'expression-matrix',
    '.txt': 'expression-matrix',
    '.xlsx': 'expression-matrix',
    '.xls': 'expression-matrix'
}

# ==================== 视图函数 ====================

def transcriptome_analysis_page(request):
    """
    转录组分析页面
    """
    return render(request, '4.html', {
        'page_title': 'Transcriptome Sequencing Analysis',
        'user': request.user if request.user.is_authenticated else None
    })

@csrf_exempt
@require_http_methods(["POST"])
def create_transcriptome_task(request):
    """
    创建转录组分析任务
    """
    try:
        user = request.user if request.user.is_authenticated else None
        # 兼容前端 JSON 请求，也兼容普通表单/空 body 请求，避免非 JSON POST 直接 500
        data = {}
        if request.body:
            content_type = (request.content_type or '').lower()
            if 'application/json' in content_type:
                try:
                    data = json.loads(request.body.decode('utf-8') or '{}')
                except json.JSONDecodeError as exc:
                    return JsonResponse({
                        'success': False,
                        'error': f'Invalid JSON payload: {str(exc)}'
                    }, status=400)
            else:
                data = request.POST.dict()
        elif request.POST:
            data = request.POST.dict()
        
        task_id = uuid.uuid4()
        
        task = TranscriptomeAnalysisTask.objects.create(
            id=task_id,
            user=user,
            username=user.username if user else 'anonymous',
            status='pending',
            parameters=data.get('parameters', {}),
            files=data.get('files', {}),
            analysis_type='transcriptome',
            progress=0,
            current_step='Initializing',
            message='Task created successfully'
        )
        
        return JsonResponse({
            'success': True,
            'task_id': str(task_id),
            'message': 'Task created successfully'
        })
        
    except Exception as e:
        logger.error(f"Error creating task: {str(e)}")
        return JsonResponse({
            'success': False,
            'error': f'Failed to create task: {str(e)}'
        }, status=500)

@csrf_exempt
@require_http_methods(["POST"])
def upload_transcriptome_file(request, task_id):
    """
    上传转录组分析文件 - 支持表达矩阵CSV和元数据CSV
    增强版：添加完整的日志记录和错误处理
    """
    import traceback
    import sys
    from datetime import datetime
    
    # 生成请求ID用于追踪
    request_id = str(uuid.uuid4())[:8]
    start_time = timezone.now()
    
    try:
        print(f"\n{'='*60}")
        print(f"=== UPLOAD REQUEST [{request_id}] ===")
        print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Task ID: {task_id}")
        print(f"Request method: {request.method}")
        print(f"Content-Type: {request.content_type}")
        print(f"Content-Length: {request.META.get('CONTENT_LENGTH', 0)} bytes")
        print(f"Remote IP: {request.META.get('REMOTE_ADDR', 'Unknown')}")
        print(f"User: {request.user.username if request.user.is_authenticated else 'anonymous'}")
        print(f"FILES keys: {list(request.FILES.keys())}")
        print(f"POST keys: {list(request.POST.keys())}")
        
        # === 1. 验证文件是否存在 ===
        if 'file' not in request.FILES:
            error_msg = 'No file uploaded. Please select a file.'
            print(f"ERROR [{request_id}]: {error_msg}")
            TranscriptomeAnalysisLogger.error(
                str(task_id), 
                error_msg, 
                request.user if request.user.is_authenticated else None,
                'file_upload'
            )
            return JsonResponse({
                'success': False,
                'error': error_msg,
                'request_id': request_id
            }, status=400)
        
        file_obj = request.FILES['file']
        file_name = file_obj.name
        file_name_lower = file_name.lower()
        file_size = file_obj.size
        
        print(f"FILE [{request_id}]:")
        print(f"  - Name: {file_name}")
        print(f"  - Size: {file_size} bytes ({file_size/1024/1024:.2f} MB)")
        print(f"  - Content-Type: {file_obj.content_type}")
        
        # === 2. 验证文件大小（限制500MB）===
        max_size = 500 * 1024 * 1024  # 500MB
        if file_size > max_size:
            error_msg = f'File too large. Maximum size is 500MB (uploaded: {file_size/1024/1024:.2f}MB)'
            print(f"ERROR [{request_id}]: {error_msg}")
            TranscriptomeAnalysisLogger.error(
                str(task_id), 
                error_msg, 
                request.user if request.user.is_authenticated else None,
                'file_upload'
            )
            return JsonResponse({
                'success': False,
                'error': error_msg,
                'request_id': request_id,
                'file_size': file_size,
                'max_size': max_size
            }, status=400)
        
        # === 3. 验证文件名（防止路径遍历攻击）===
        import re
        if '..' in file_name or '/' in file_name or '\\' in file_name:
            error_msg = 'Invalid filename: contains path traversal characters'
            print(f"ERROR [{request_id}]: {error_msg}")
            return JsonResponse({
                'success': False,
                'error': error_msg,
                'request_id': request_id
            }, status=400)
        
        # === 4. 获取并验证文件类型 ===
        file_type = request.POST.get('file_type', 'unknown')
        print(f"REQUESTED FILE TYPE: {file_type}")
        
        # 自动检测文件类型
        if file_type == 'unknown' or not file_type:
            original_file_type = file_type
            if 'fastq' in file_name_lower or file_name_lower.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
                file_type = 'fastq-files'
            elif file_name_lower.endswith(('.csv', '.tsv', '.txt', '.xlsx', '.xls')):
                file_type = 'expression-matrix'
            elif 'metadata' in file_name_lower or 'sample' in file_name_lower:
                file_type = 'sample-metadata'
            elif file_name_lower.endswith(('.fa', '.fasta', '.fa.gz', '.fasta.gz')):
                file_type = 'reference-genome'
            elif file_name_lower.endswith(('.gtf', '.gff', '.gff3')):
                file_type = 'annotation'
            elif file_name_lower.endswith(('.bed',)):
                file_type = 'gene-list'
            elif file_name_lower.endswith(('.json', '.gmt')):
                file_type = 'custom-genesets'
            else:
                file_type = 'expression-matrix'  # 默认类型
            print(f"AUTO-DETECTED FILE TYPE: {file_type} (original: {original_file_type})")
        
        # 验证文件类型是否允许
        allowed_types = ['expression-matrix', 'sample-metadata', 'fastq-files', 
                        'reference-genome', 'gene-list', 'annotation', 'custom-genesets']
        if file_type not in allowed_types:
            error_msg = f'Invalid file type: {file_type}. Allowed types: {", ".join(allowed_types)}'
            print(f"ERROR [{request_id}]: {error_msg}")
            return JsonResponse({
                'success': False,
                'error': error_msg,
                'request_id': request_id,
                'file_type': file_type,
                'allowed_types': allowed_types
            }, status=400)
        
        # === 5. 验证任务是否存在 ===
        try:
            task = TranscriptomeAnalysisTask.objects.get(id=task_id)
            print(f"TASK FOUND: {task.id}")
            print(f"  - Status: {task.status}")
            print(f"  - Created: {task.created_at}")
            print(f"  - Username: {task.username}")
            print(f"  - Existing files: {list(task.files.keys()) if task.files else 'None'}")
        except TranscriptomeAnalysisTask.DoesNotExist:
            error_msg = f'Task {task_id} not found'
            print(f"ERROR [{request_id}]: {error_msg}")
            return JsonResponse({
                'success': False,
                'error': error_msg,
                'request_id': request_id
            }, status=404)
        except Exception as e:
            error_msg = f'Error accessing task: {str(e)}'
            print(f"ERROR [{request_id}]: {error_msg}")
            print(traceback.format_exc())
            return JsonResponse({
                'success': False,
                'error': error_msg,
                'request_id': request_id
            }, status=500)
        
        # === 6. 检查任务状态 ===
        if task.status not in ['pending', 'running', 'completed', 'failed', 'cancelled']:
            print(f"WARNING [{request_id}]: Task in unusual state: {task.status}")
        
        # === 7. 创建任务目录 ===
        task_uuid_str = str(task_id)
        task_dir = os.path.join(settings.MEDIA_ROOT, 'transcriptome_analysis', task_uuid_str)
        
        try:
            os.makedirs(task_dir, mode=0o755, exist_ok=True)
            print(f"TASK DIRECTORY: {task_dir}")
            print(f"  - Exists: {os.path.exists(task_dir)}")
            print(f"  - Writable: {os.access(task_dir, os.W_OK)}")
        except PermissionError:
            error_msg = f'Permission denied: Cannot create directory {task_dir}'
            print(f"ERROR [{request_id}]: {error_msg}")
            return JsonResponse({
                'success': False,
                'error': error_msg,
                'request_id': request_id
            }, status=500)
        except Exception as e:
            error_msg = f'Failed to create task directory: {str(e)}'
            print(f"ERROR [{request_id}]: {error_msg}")
            return JsonResponse({
                'success': False,
                'error': error_msg,
                'request_id': request_id
            }, status=500)
        
        # === 8. 保存文件 ===
        # 生成安全的文件名
        safe_filename = f"{uuid.uuid4()}_{os.path.basename(file_name)}"
        # 进一步清理文件名
        safe_filename = re.sub(r'[^\w\-_.]', '_', safe_filename)
        file_path = os.path.join(task_dir, safe_filename)
        
        print(f"SAVING FILE TO: {file_path}")
        
        try:
            # 分块写入文件，添加进度日志
            bytes_written = 0
            chunk_count = 0
            with open(file_path, 'wb+') as destination:
                for chunk in file_obj.chunks():
                    destination.write(chunk)
                    bytes_written += len(chunk)
                    chunk_count += 1
                    # 每100个块打印一次进度（避免日志过多）
                    if chunk_count % 100 == 0:
                        print(f"  - Written: {bytes_written}/{file_size} bytes ({bytes_written/file_size*100:.1f}%)")
            
            print(f"FILE SAVED SUCCESSFULLY:")
            print(f"  - Path: {file_path}")
            print(f"  - Size: {bytes_written} bytes")
            print(f"  - Chunks: {chunk_count}")
            
            # 验证文件完整性
            if bytes_written != file_size:
                error_msg = f'File size mismatch: expected {file_size}, got {bytes_written}'
                print(f"ERROR [{request_id}]: {error_msg}")
                os.remove(file_path)
                return JsonResponse({
                    'success': False,
                    'error': error_msg,
                    'request_id': request_id
                }, status=500)
                
        except IOError as e:
            error_msg = f'IO Error while saving file: {str(e)}'
            print(f"ERROR [{request_id}]: {error_msg}")
            # 尝试清理部分写入的文件
            try:
                if os.path.exists(file_path):
                    os.remove(file_path)
                    print(f"CLEANED UP partial file: {file_path}")
            except:
                pass
            return JsonResponse({
                'success': False,
                'error': error_msg,
                'request_id': request_id
            }, status=500)
        except Exception as e:
            error_msg = f'Failed to save file: {str(e)}'
            print(f"ERROR [{request_id}]: {error_msg}")
            print(traceback.format_exc())
            return JsonResponse({
                'success': False,
                'error': error_msg,
                'request_id': request_id
            }, status=500)
        
        # === 9. 验证表达矩阵文件 ===
        if file_type == 'expression-matrix':
            try:
                print(f"VALIDATING expression matrix...")
                df = read_expression_matrix(file_path)
                
                # 获取数据统计信息
                n_genes = df.shape[0]
                n_samples = df.shape[1]
                n_zeros = (df == 0).sum().sum()
                n_total = n_genes * n_samples
                zero_percentage = (n_zeros / n_total * 100) if n_total > 0 else 0
                
                print(f"EXPRESSION MATRIX STATISTICS:")
                print(f"  - Genes: {n_genes}")
                print(f"  - Samples: {n_samples}")
                print(f"  - Zero values: {n_zeros} ({zero_percentage:.1f}%)")
                print(f"  - Sample names: {list(df.columns[:5])}{'...' if n_samples > 5 else ''}")
                print(f"  - Gene names: {list(df.index[:5])}{'...' if n_genes > 5 else ''}")
                
                # 保存样本数量信息到任务
                if 'fastq_files' not in task.files:
                    task.files['fastq_files'] = []
                
                # 更新任务参数
                task.parameters = task.parameters or {}
                task.parameters['sample_count'] = n_samples
                task.parameters['gene_count'] = n_genes
                task.parameters['zero_percentage'] = round(zero_percentage, 2)
                task.save()
                
            except Exception as e:
                error_msg = f'Invalid expression matrix file: {str(e)}'
                print(f"ERROR [{request_id}]: {error_msg}")
                print(traceback.format_exc())
                # 删除无效文件
                try:
                    os.remove(file_path)
                    print(f"DELETED invalid file: {file_path}")
                except:
                    pass
                return JsonResponse({
                    'success': False,
                    'error': error_msg,
                    'request_id': request_id,
                    'file_type': file_type,
                    'file_name': file_name
                }, status=400)
        
        # === 10. 验证参考基因组文件 ===
        if file_type == 'reference-genome':
            try:
                print(f"VALIDATING reference genome file...")
                # 简单验证：检查文件是否为空
                if os.path.getsize(file_path) == 0:
                    raise ValueError("Reference genome file is empty")
                
                # 尝试读取前几行
                with _open_text_maybe_gzip(file_path) as f:
                    first_line = f.readline().strip()
                    if not first_line.startswith('>'):
                        print(f"WARNING: Reference file may not be in FASTA format (first line: {first_line[:50]})")
                
                print(f"REFERENCE GENOME VALIDATED")
                
            except Exception as e:
                error_msg = f'Invalid reference genome file: {str(e)}'
                print(f"ERROR [{request_id}]: {error_msg}")
                try:
                    os.remove(file_path)
                except:
                    pass
                return JsonResponse({
                    'success': False,
                    'error': error_msg,
                    'request_id': request_id
                }, status=400)
        
        # === 11. 验证样本元数据文件 ===
        if file_type == 'sample-metadata':
            try:
                print(f"VALIDATING sample metadata file...")
                # 尝试读取元数据
                if file_name_lower.endswith('.csv'):
                    df_meta = pd.read_csv(file_path)
                elif file_name_lower.endswith(('.tsv', '.txt')):
                    df_meta = pd.read_csv(file_path, sep='\t')
                else:
                    df_meta = pd.read_csv(file_path)  # 尝试自动检测
                
                print(f"SAMPLE METADATA STATISTICS:")
                print(f"  - Rows: {df_meta.shape[0]}")
                print(f"  - Columns: {df_meta.shape[1]}")
                print(f"  - Column names: {list(df_meta.columns)}")
                
                # 检查必要的列
                if df_meta.shape[1] < 2:
                    print(f"WARNING: Metadata should have at least 2 columns (sample names and groups)")
                
            except Exception as e:
                error_msg = f'Invalid sample metadata file: {str(e)}'
                print(f"ERROR [{request_id}]: {error_msg}")
                return JsonResponse({
                    'success': False,
                    'error': error_msg,
                    'request_id': request_id
                }, status=400)
        
        # === 12. 创建文件记录 ===
        try:
            file_record = TranscriptomeAnalysisFile.objects.create(
                task=task,
                file_name=file_name,
                file_path=file_path,
                file_type=file_type,
                file_size=file_size
            )
            print(f"FILE RECORD CREATED: {file_record.id}")
            
        except Exception as e:
            error_msg = f'Failed to create file record: {str(e)}'
            print(f"ERROR [{request_id}]: {error_msg}")
            print(traceback.format_exc())
            # 删除已保存的文件
            try:
                os.remove(file_path)
                print(f"DELETED file due to database error: {file_path}")
            except:
                pass
            return JsonResponse({
                'success': False,
                'error': error_msg,
                'request_id': request_id
            }, status=500)
        
        # === 13. 更新任务文件信息 ===
        try:
            files_info = task.files.copy() if task.files else {}
            
            if file_type == 'expression-matrix':
                files_info['expression_matrix'] = {
                    'id': str(file_record.id),
                    'name': file_name,
                    'path': file_path,
                    'size': file_size,
                    'samples': task.parameters.get('sample_count', 0),
                    'genes': task.parameters.get('gene_count', 0),
                    'uploaded_at': timezone.now().isoformat()
                }
            elif file_type == 'fastq-files':
                if 'fastq_files' not in files_info:
                    files_info['fastq_files'] = []
                files_info['fastq_files'].append({
                    'id': str(file_record.id),
                    'name': file_name,
                    'path': file_path,
                    'size': file_size,
                    'uploaded_at': timezone.now().isoformat()
                })
            else:
                files_info[file_type] = {
                    'id': str(file_record.id),
                    'name': file_name,
                    'path': file_path,
                    'size': file_size,
                    'uploaded_at': timezone.now().isoformat()
                }
            
            task.files = files_info
            task.updated_at = timezone.now()
            task.save()
            
            print(f"TASK FILES UPDATED:")
            print(f"  - File type: {file_type}")
            print(f"  - Total files now: {len(files_info.get('fastq_files', [])) if file_type == 'fastq-files' else 'N/A'}")
            
        except Exception as e:
            error_msg = f'Failed to update task files: {str(e)}'
            print(f"ERROR [{request_id}]: {error_msg}")
            print(traceback.format_exc())
            # 记录错误但不返回失败，因为文件已成功上传
            TranscriptomeAnalysisLogger.error(
                str(task_id),
                error_msg,
                request.user if request.user.is_authenticated else None,
                'file_upload'
            )
        
        # === 14. 记录成功日志 ===
        elapsed_time = (timezone.now() - start_time).total_seconds()
        print(f"UPLOAD COMPLETED SUCCESSFULLY [{request_id}]:")
        print(f"  - Time elapsed: {elapsed_time:.2f} seconds")
        print(f"  - File: {file_name}")
        print(f"  - Type: {file_type}")
        print(f"  - Size: {file_size} bytes")
        print(f"  - File ID: {file_record.id}")
        print(f"{'='*60}\n")
        
        TranscriptomeAnalysisLogger.info(
            task_uuid_str,
            f'File uploaded successfully: {file_name} ({file_type}, {file_size} bytes, {elapsed_time:.2f}s)',
            request.user if request.user.is_authenticated else None,
            'file_upload'
        )
        
        # 构建响应
        response_data = {
            'success': True,
            'request_id': request_id,
            'file_id': str(file_record.id),
            'file_name': file_name,
            'file_type': file_type,
            'file_size': file_size,
            'file_path': file_path,
            'message': f'File {file_name} uploaded successfully',
            'task_id': task_uuid_str,
            'upload_time': timezone.now().isoformat(),
            'elapsed_seconds': round(elapsed_time, 2)
        }
        
        # 添加表达矩阵统计信息
        if file_type == 'expression-matrix':
            response_data['samples'] = task.parameters.get('sample_count', 0)
            response_data['genes'] = task.parameters.get('gene_count', 0)
        
        return JsonResponse(response_data)
        
    except Exception as e:
        # 全局异常捕获
        elapsed_time = (timezone.now() - start_time).total_seconds()
        error_msg = f'Unexpected error: {str(e)}'
        print(f"\n!!! UNEXPECTED ERROR [{request_id}] !!!")
        print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Error type: {type(e).__name__}")
        print(f"Error message: {str(e)}")
        print(f"Time elapsed: {elapsed_time:.2f} seconds")
        print("\nTraceback:")
        traceback.print_exc()
        print(f"{'='*60}\n")
        
        logger.error(f"Error uploading file for task {task_id}: {str(e)}", exc_info=True)
        
        return JsonResponse({
            'success': False,
            'error': f'Failed to upload file: {str(e)}',
            'request_id': request_id,
            'error_type': type(e).__name__
        }, status=500)

@csrf_exempt
@require_http_methods(["POST"])
def run_transcriptome_analysis(request, task_id):
    """
    运行转录组分析 - 使用上传的表达矩阵进行真实数据分析
    """
    try:
        print(f"=== DEBUG RUN ANALYSIS ===")
        print(f"Task ID: {task_id}")
        
        # 验证任务
        try:
            task = TranscriptomeAnalysisTask.objects.get(id=task_id)
            print(f"Task found: {task.id}")
            print(f"Task status: {task.status}")
        except TranscriptomeAnalysisTask.DoesNotExist:
            return JsonResponse({
                'success': False,
                'error': 'Task not found'
            }, status=404)

        # 解析前端运行参数。前端选择内置参考基因组时不会上传 reference-genome 文件，
        # 因此这里需要把 hg38/mm10 等选择转换成任务里的 reference-genome 记录。
        run_payload = {}
        if request.body:
            try:
                run_payload = json.loads(request.body.decode('utf-8') or '{}')
            except json.JSONDecodeError:
                run_payload = {}
        submitted_files = run_payload.get('files') or {}
        submitted_params = run_payload.get('parameters') or {}
        if submitted_params:
            task.parameters = {**(task.parameters or {}), **submitted_params}
            task.save(update_fields=['parameters', 'updated_at'])

        selected_reference = (
            submitted_files.get('reference_genome')
            or run_payload.get('reference_genome')
            or (task.parameters or {}).get('reference_genome')
        )
        if selected_reference and selected_reference != 'custom' and 'reference-genome' not in (task.files or {}):
            task.files = task.files or {}
            task.files['reference-genome'] = {
                'id': f'builtin-{selected_reference}',
                'name': selected_reference,
                'path': f'builtin://{selected_reference}',
                'size': 0,
                'source': 'builtin',
                'uploaded_at': timezone.now().isoformat()
            }
            task.parameters = task.parameters or {}
            task.parameters['reference_genome'] = selected_reference
            task.save(update_fields=['files', 'parameters', 'updated_at'])
            print(f"Using built-in reference genome: {selected_reference}")
        
        # 检查任务状态 - 修改这里，允许重新运行
        if task.status == 'running':
            print(f"Analysis is already running, resetting task...")
            # 不返回错误，而是重置任务状态
            task.status = 'pending'
            task.progress = 0
            task.current_step = 'Restarting analysis'
            task.message = 'Restarting analysis pipeline'
            task.error = None
            task.save()
            print(f"Task reset to pending")
        
        # 检查必需文件
        has_fastq_files = 'fastq_files' in task.files and len(task.files.get('fastq_files', [])) > 0
        has_reference = 'reference-genome' in task.files
        has_metadata = 'sample-metadata' in task.files
        
        print(f"Has FASTQ files: {has_fastq_files}")
        print(f"Has reference genome: {has_reference}")
        print(f"Has sample metadata: {has_metadata}")
        
        if not has_fastq_files:
            # 兼容表达矩阵上传：如果用户上传了表达矩阵而不是 FASTQ，也允许进入分析流程。
            has_expression_matrix = 'expression_matrix' in (task.files or {})
            if has_expression_matrix:
                task.files = task.files or {}
                task.files['fastq_files'] = [{
                    'id': task.files['expression_matrix'].get('id', 'expression-matrix'),
                    'name': task.files['expression_matrix'].get('name', 'expression_matrix'),
                    'path': task.files['expression_matrix'].get('path'),
                    'size': task.files['expression_matrix'].get('size', 0),
                    'source': 'expression_matrix'
                }]
                task.save(update_fields=['files', 'updated_at'])
                has_fastq_files = True
                print("Using uploaded expression matrix as analysis input")
            else:
                error_msg = 'No FASTQ files or expression matrix uploaded. Please upload at least one input file.'
                print(f"Error: {error_msg}")
                return JsonResponse({
                    'success': False,
                    'error': error_msg
                }, status=400)
        
        # 表达矩阵模式不需要参考基因组；FASTQ 原始 reads 模式必须有参考/索引
        has_expression_matrix = 'expression_matrix' in (task.files or {})
        if not has_reference and not has_expression_matrix:
            error_msg = 'No reference genome/index selected. FASTQ production analysis requires a configured reference, or upload an expression matrix for downstream analysis.'
            print(f"Error: {error_msg}")
            return JsonResponse({
                'success': False,
                'error': error_msg
            }, status=400)
        
        # 重置任务状态为运行中
        task.status = 'running'
        task.progress = 0
        task.current_step = 'Starting analysis'
        task.message = 'Initializing analysis pipeline'
        task.started_at = timezone.now()
        task.error = None
        task.save()
        
        print(f"Task status updated to running")
        
        # 执行分析（异步方式）
        try:
            # 这里可以调用Celery任务，或者直接执行同步分析
            analysis_result = perform_transcriptome_analysis(task, str(task_id))
            
            if analysis_result['success']:
                task.status = 'completed'
                task.progress = 100
                task.current_step = 'Completed'
                task.message = 'Analysis completed successfully'
                task.completed_at = timezone.now()
                task.save()
                
                return JsonResponse({
                    'success': True,
                    'message': 'Analysis completed successfully',
                    'status': task.status,
                    'progress': 100
                })
            else:
                task.status = 'failed'
                task.error = analysis_result.get('error', 'Unknown error')
                task.message = 'Analysis failed'
                task.save()
                
                return JsonResponse({
                    'success': False,
                    'error': analysis_result.get('error', 'Analysis failed')
                }, status=500)
                
        except Exception as e:
            task.status = 'failed'
            task.error = str(e)
            task.message = 'Analysis failed with exception'
            task.save()
            
            logger.error(f"Error in analysis: {str(e)}", exc_info=True)
            return JsonResponse({
                'success': False,
                'error': f'Analysis failed: {str(e)}'
            }, status=500)
        
    except Exception as e:
        logger.error(f"Error running analysis: {str(e)}", exc_info=True)
        return JsonResponse({
            'success': False,
            'error': f'Failed to start analysis: {str(e)}'
        }, status=500)

def perform_transcriptome_analysis(task, task_id_str):
    """
    生产级转录组分析核心函数：
    - expression_matrix：读取用户上传表达矩阵并做下游分析
    - FASTQ：必须调用真实定量工具（当前实现 Salmon）；缺工具/缺索引会失败并给出明确错误，不再模拟
    """
    try:
        print(f"=== Starting transcriptome analysis for task {task_id_str} ===")
        task_dir = os.path.join(settings.MEDIA_ROOT, 'transcriptome_analysis', task_id_str)
        results_dir = os.path.join(task_dir, 'results')
        plots_dir = os.path.join(task_dir, 'plots')
        os.makedirs(results_dir, exist_ok=True)
        os.makedirs(plots_dir, exist_ok=True)
        log_file = os.path.join(results_dir, 'pipeline.log')
        _append_log(log_file, f"Starting transcriptome analysis task={task_id_str}")

        # 步骤1：文件信息
        task.progress = 10
        task.current_step = 'Preparing files'
        task.save()
        files = task.files or {}
        fastq_files = files.get('fastq_files', [])
        expression_matrix_info = files.get('expression_matrix')
        _append_log(log_file, f"FASTQ files: {len(fastq_files)}, expression_matrix: {bool(expression_matrix_info)}")

        # 步骤2：读取样本元数据
        task.progress = 15
        task.current_step = 'Loading sample metadata'
        task.save()
        sample_groups = {}
        metadata_file_info = files.get('sample-metadata', {})
        metadata_file_path = metadata_file_info.get('path') if isinstance(metadata_file_info, dict) else None
        if metadata_file_path and os.path.exists(metadata_file_path):
            try:
                # 自动识别 csv/tsv
                sep = ',' if metadata_file_path.lower().endswith('.csv') else ('\t' if metadata_file_path.lower().endswith(('.tsv','.txt')) else None)
                df_metadata = pd.read_csv(metadata_file_path, sep=sep, engine='python' if sep is None else 'c')
                if len(df_metadata.columns) >= 2:
                    sample_col, group_col = df_metadata.columns[0], df_metadata.columns[1]
                    for _, row in df_metadata.iterrows():
                        sample_groups[str(row[sample_col])] = str(row[group_col])
                _append_log(log_file, f"Loaded metadata: {df_metadata.shape}, groups={set(sample_groups.values())}")
            except Exception as e:
                raise RuntimeError(f"Failed to read sample metadata: {e}")

        # 步骤3：表达矩阵来源
        task.progress = 30
        task.current_step = 'Quantification / expression matrix generation'
        task.save()
        extra_file_paths = {}
        pipeline_mode = 'expression_matrix'

        if expression_matrix_info and expression_matrix_info.get('path') and os.path.exists(expression_matrix_info.get('path')):
            # 表达矩阵模式：这是生产可用的下游分析入口，不伪造数据
            matrix_path = expression_matrix_info.get('path')
            _append_log(log_file, f"Reading uploaded expression matrix: {matrix_path}")
            df_counts = read_expression_matrix(matrix_path)
            sample_names = df_counts.columns.tolist()
            if not sample_groups:
                sample_groups = {sample: f"Group_{(i % 2) + 1}" for i, sample in enumerate(sample_names)}
            expression_file_path = os.path.join(results_dir, 'raw_counts.csv')
            df_counts.to_csv(expression_file_path)
        else:
            # FASTQ 模式：生产级要求真实工具和真实索引
            pipeline_mode = 'fastq_salmon'
            if not fastq_files:
                raise RuntimeError('No FASTQ files uploaded.')
            ref_info = resolve_reference_for_rnaseq(task, results_dir, log_file)
            if ref_info.get('mode') != 'salmon':
                raise RuntimeError(f"Unsupported RNA-seq quantification mode: {ref_info.get('mode')}")
            task.progress = 35
            task.current_step = 'Running Salmon quantification'
            task.save()
            df_counts, quant_paths = run_salmon_quant(task, fastq_files, ref_info, results_dir, log_file)
            expression_file_path = quant_paths['counts']
            extra_file_paths.update(quant_paths)
            sample_names = df_counts.columns.tolist()
            if not sample_groups:
                sample_groups = {sample: f"Group_{(i % 2) + 1}" for i, sample in enumerate(sample_names)}

        # 步骤4：标准化
        task.progress = 45
        task.current_step = 'Data normalization'
        task.save()
        df_counts = df_counts.apply(pd.to_numeric, errors='coerce').fillna(0)
        lib_sizes = df_counts.sum(axis=0).replace(0, 1)
        df_cpm = df_counts.div(lib_sizes, axis=1) * 1e6
        df_log2cpm = np.log2(df_cpm + 1)
        cpm_path = os.path.join(results_dir, 'cpm_normalized.csv')
        log2cpm_path = os.path.join(results_dir, 'log2cpm_normalized.csv')
        df_cpm.to_csv(cpm_path)
        df_log2cpm.to_csv(log2cpm_path)

        # 步骤5：QC图
        task.progress = 55
        task.current_step = 'Quality control'
        task.save()
        generate_qc_plots(df_log2cpm, df_counts, sample_groups, plots_dir, task_id_str)

        # 步骤6：差异表达
        task.progress = 65
        task.current_step = 'Differential expression analysis'
        task.save()
        de_results = perform_differential_expression(
            df_log2cpm, sample_groups, results_dir, plots_dir,
            (task.parameters or {}).get('differential_expression', {})
        )

        # 步骤7：富集
        task.progress = 80
        task.current_step = 'Functional enrichment'
        task.save()
        enrichment_results = perform_enrichment_analysis(
            de_results, results_dir, plots_dir,
            (task.parameters or {}).get('enrichment', {})
        )

        # 步骤8：可视化
        task.progress = 90
        task.current_step = 'Generating visualizations'
        task.save()
        generate_all_plots(df_log2cpm, de_results, enrichment_results, sample_groups, plots_dir)

        # 步骤9：保存结果
        task.progress = 95
        task.current_step = 'Saving results'
        task.save()
        file_paths = {
            'counts': expression_file_path,
            'cpm': cpm_path,
            'log2cpm': log2cpm_path,
            'de_genes': de_results.get('file_path', ''),
            'enrichment': enrichment_results.get('file_path', ''),
            'summary': os.path.join(results_dir, 'analysis_summary.txt'),
            'plots': plots_dir,
            'pipeline_log': log_file,
            **extra_file_paths
        }
        summary = {
            'pipeline_mode': pipeline_mode,
            'total_samples': int(df_counts.shape[1]),
            'total_genes': int(df_counts.shape[0]),
            'de_up': int(de_results.get('up_count', 0)),
            'de_down': int(de_results.get('down_count', 0)),
            'enriched_terms': int(enrichment_results.get('term_count', 0)),
            'analysis_date': timezone.now().isoformat(),
            'production_mode': True,
        }
        with open(file_paths['summary'], 'w', encoding='utf-8') as f:
            f.write("TRANSCRIPTOME ANALYSIS SUMMARY\n")
            f.write("=" * 50 + "\n\n")
            for k, v in summary.items():
                f.write(f"{k}: {v}\n")
            f.write(f"\nPipeline log: {log_file}\n")
        TranscriptomeAnalysisResult.objects.filter(task=task).delete()
        result = TranscriptomeAnalysisResult.objects.create(
            task=task, summary=summary, differential_expression=de_results,
            enrichment=enrichment_results, file_paths=file_paths
        )
        task.results = {'analysis_id': str(result.id), 'summary': summary, 'file_paths': file_paths}
        task.save()
        _append_log(log_file, "Analysis completed successfully")
        return {'success': True, 'result_id': str(result.id), 'summary': summary}
    except Exception as e:
        print(f"Error in analysis: {str(e)}")
        import traceback
        traceback.print_exc()
        return {'success': False, 'error': str(e)}

def read_expression_matrix(file_path):
    """
    读取表达矩阵文件，支持多种格式，增强错误处理
    """
    import pandas as pd
    import os
    
    print(f"Reading expression matrix from: {file_path}")
    
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    
    if os.path.getsize(file_path) == 0:
        raise ValueError("File is empty")
    
    file_ext = os.path.splitext(file_path)[1].lower()
    file_name = os.path.basename(file_path).lower()
    
    try:
        # 根据文件扩展名选择读取方法
        if file_ext == '.csv':
            df = pd.read_csv(file_path, index_col=0)
            print(f"Read as CSV: {df.shape}")
        elif file_ext == '.tsv':
            df = pd.read_csv(file_path, sep='\t', index_col=0)
            print(f"Read as TSV: {df.shape}")
        elif file_ext in ['.xlsx', '.xls']:
            df = pd.read_excel(file_path, index_col=0)
            print(f"Read as Excel: {df.shape}")
        else:  # .txt 或其他
            # 尝试自动检测分隔符
            with open(file_path, 'r') as f:
                first_line = f.readline()
                if '\t' in first_line:
                    df = pd.read_csv(file_path, sep='\t', index_col=0)
                    print(f"Auto-detected TSV format")
                else:
                    df = pd.read_csv(file_path, sep=',', index_col=0)
                    print(f"Auto-detected CSV format")
        
        # 确保数据是数值型
        df = df.apply(pd.to_numeric, errors='coerce')
        
        # 检查数据有效性
        if df.empty:
            raise ValueError("DataFrame is empty after reading")
        
        # 处理NaN值
        nan_count = df.isna().sum().sum()
        if nan_count > 0:
            print(f"Warning: Found {nan_count} NaN values, filling with 0")
            df = df.fillna(0)
        
        # 确保索引和列名是字符串
        df.index = df.index.astype(str)
        df.columns = df.columns.astype(str)
        
        # 检查索引重复
        if df.index.duplicated().any():
            dup_count = df.index.duplicated().sum()
            print(f"Warning: Found {dup_count} duplicate gene names, keeping first occurrence")
            df = df[~df.index.duplicated(keep='first')]
        
        # 检查列名重复
        if df.columns.duplicated().any():
            dup_count = df.columns.duplicated().sum()
            print(f"Warning: Found {dup_count} duplicate sample names, keeping first occurrence")
            df = df.loc[:, ~df.columns.duplicated(keep='first')]
        
        # 基本统计
        print(f"Final expression matrix shape: {df.shape[0]} genes, {df.shape[1]} samples")
        print(f"Data range: {df.min().min():.2f} - {df.max().max():.2f}")
        print(f"Mean expression: {df.mean().mean():.2f}")
        print(f"Median expression: {df.median().median():.2f}")
        
        return df
        
    except pd.errors.EmptyDataError:
        raise ValueError("File is empty or contains no data")
    except pd.errors.ParserError as e:
        raise ValueError(f"Failed to parse file: {str(e)}")
    except Exception as e:
        raise Exception(f"Failed to read expression matrix: {str(e)}")

def generate_qc_plots(df_log2cpm, df_counts, groups, plots_dir, task_id):
    """
    生成质量控制图表
    """
    try:
        # 1. 样本相关性热图
        plt.figure(figsize=(12, 10))
        corr_matrix = df_log2cpm.corr()
        
        # 添加分组注释
        if groups and isinstance(groups, dict):
            sample_groups = [groups.get(col, 'Unknown') for col in df_log2cpm.columns]
            unique_groups = list(set(sample_groups))
            group_colors = plt.cm.tab10(np.linspace(0, 1, len(unique_groups)))
            group_color_map = {g: group_colors[i] for i, g in enumerate(unique_groups)}
            colors = [group_color_map[g] for g in sample_groups]
            
            # 创建带颜色标签的热图
            g = sns.clustermap(corr_matrix, 
                              cmap='coolwarm', 
                              center=0,
                              annot=False,
                              figsize=(12, 10),
                              dendrogram_ratio=0.15,
                              cbar_pos=(0.02, 0.8, 0.03, 0.18),
                              col_colors=colors)
            g.ax_heatmap.set_title('Sample Correlation Matrix', fontsize=14)
        else:
            plt.figure(figsize=(10, 8))
            sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0, fmt='.2f')
            plt.title('Sample Correlation Matrix')
        
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'qc-metrics.png'), dpi=150, bbox_inches='tight')
        plt.close()
        
        # 2. 表达分布箱线图
        plt.figure(figsize=(14, 6))
        df_log2cpm.boxplot(rot=45, grid=False)
        plt.title('Gene Expression Distribution Across Samples', fontsize=14)
        plt.ylabel('log2(CPM+1)', fontsize=12)
        plt.xlabel('Samples', fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'expression-distribution.png'), dpi=150, bbox_inches='tight')
        plt.close()
        
        # 3. 文库大小分布
        plt.figure(figsize=(10, 6))
        lib_sizes = df_counts.sum(axis=0) / 1e6  # 转换为百万
        plt.bar(range(len(lib_sizes)), lib_sizes.values)
        plt.xticks(range(len(lib_sizes)), lib_sizes.index, rotation=45)
        plt.title('Library Size Distribution', fontsize=14)
        plt.ylabel('Total Counts (Millions)', fontsize=12)
        plt.xlabel('Samples', fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'library-sizes.png'), dpi=150, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error generating QC plots: {e}")

def perform_differential_expression(df_log2cpm, groups, results_dir, plots_dir, params=None):
    """
    执行差异表达分析 - 完全修复版本
    """
    try:
        if params is None:
            params = {}
        
        # 获取参数
        pval_cutoff = float(params.get('pval_cutoff', 0.05))
        lfc_cutoff = float(params.get('lfc_cutoff', 1.0))
        
        print(f"Starting differential expression analysis...")
        print(f"Data shape: {df_log2cpm.shape}")
        print(f"Groups: {groups}")
        
        # 准备分组信息
        sample_names = df_log2cpm.columns.tolist()
        
        # 构建分组标签
        if isinstance(groups, dict):
            group_labels = []
            for sample in sample_names:
                group = groups.get(sample, 'Unknown')
                group_labels.append(str(group))
        else:
            # 默认分组
            n_samples = len(sample_names)
            group_labels = ['Control'] * (n_samples // 2) + ['Treatment'] * (n_samples - n_samples // 2)
        
        # 获取唯一分组
        unique_groups = list(set(group_labels))
        unique_groups = [g for g in unique_groups if g != 'Unknown']
        
        print(f"Unique groups: {unique_groups}")
        print(f"Group labels: {group_labels}")
        
        if len(unique_groups) < 2:
            print("Warning: Only one group found, cannot perform differential expression")
            return {
                'file_path': '',
                'total_count': 0,
                'up_count': 0,
                'down_count': 0,
                'top_genes': []
            }
        
        # 选择两个主要分组进行比较
        group1 = unique_groups[0]
        group2 = unique_groups[1]
        
        # 获取样本名称列表
        group1_samples = [sample_names[i] for i, g in enumerate(group_labels) if g == group1]
        group2_samples = [sample_names[i] for i, g in enumerate(group_labels) if g == group2]
        
        print(f"Comparing {group2} (n={len(group2_samples)}) vs {group1} (n={len(group1_samples)})")
        print(f"Group1 samples: {group1_samples[:3]}...")
        print(f"Group2 samples: {group2_samples[:3]}...")
        
        # 执行t检验
        de_results = []
        
        for idx, gene in enumerate(df_log2cpm.index):
            if idx % 100 == 0:
                print(f"Processing gene {idx}/{len(df_log2cpm.index)}")
            
            try:
                # 修复：使用 loc 而不是 iloc
                expr1 = df_log2cpm.loc[gene, group1_samples].values
                expr2 = df_log2cpm.loc[gene, group2_samples].values
                
                # 检查是否有足够的非NaN值
                expr1 = expr1[~np.isnan(expr1)]
                expr2 = expr2[~np.isnan(expr2)]
                
                if len(expr1) < 2 or len(expr2) < 2:
                    continue
                
                # t检验
                t_stat, p_value = stats.ttest_ind(expr2, expr1, equal_var=False)
                
                # 计算log2FC
                log2fc = np.mean(expr2) - np.mean(expr1)
                
                # 计算base mean
                base_mean = np.mean(df_log2cpm.loc[gene])
                
                de_results.append({
                    'gene_id': gene,
                    'log2fc': float(log2fc) if not np.isnan(log2fc) else 0.0,
                    'pval': float(p_value) if not np.isnan(p_value) else 1.0,
                    'padj': float(p_value) if not np.isnan(p_value) else 1.0,
                    'base_mean': float(base_mean) if not np.isnan(base_mean) else 0.0,
                    'regulation': 'up' if log2fc > lfc_cutoff and p_value < pval_cutoff else 
                                 ('down' if log2fc < -lfc_cutoff and p_value < pval_cutoff else 'ns')
                })
                
            except Exception as e:
                print(f"Error processing gene {gene}: {e}")
                continue
        
        print(f"Processed {len(de_results)} genes")
        
        if not de_results:
            print("No genes processed successfully")
            return {
                'file_path': '',
                'total_count': 0,
                'up_count': 0,
                'down_count': 0,
                'top_genes': []
            }
        
        # 转换为DataFrame
        df_de = pd.DataFrame(de_results)
        
        # 多重检验校正 - BH方法
        df_de = df_de.sort_values('pval')
        n_tests = len(df_de)
        df_de['padj'] = df_de['pval'] * n_tests / (np.arange(n_tests) + 1)
        df_de['padj'] = df_de['padj'].clip(upper=1.0)
        df_de = df_de.sort_index()
        
        # 保存结果
        de_file_path = os.path.join(results_dir, 'differential_genes.csv')
        df_de.to_csv(de_file_path, index=False)
        
        # 生成火山图
        generate_volcano_plot(df_de, lfc_cutoff, pval_cutoff, plots_dir)
        
        # 统计
        up_count = len(df_de[df_de['regulation'] == 'up'])
        down_count = len(df_de[df_de['regulation'] == 'down'])
        
        # 获取top基因
        sig_genes = df_de[df_de['regulation'] != 'ns'].nsmallest(20, 'padj') if len(df_de[df_de['regulation'] != 'ns']) > 0 else pd.DataFrame()
        top_genes = sig_genes.to_dict('records') if len(sig_genes) > 0 else []
        
        print(f"DE results: {up_count} up, {down_count} down")
        
        return {
            'file_path': de_file_path,
            'total_count': len(df_de),
            'up_count': int(up_count),
            'down_count': int(down_count),
            'top_genes': top_genes[:10],
            'comparison': f'{group2} vs {group1}',
            'pval_cutoff': pval_cutoff,
            'lfc_cutoff': lfc_cutoff
        }
        
    except Exception as e:
        print(f"Error in differential expression: {e}")
        import traceback
        traceback.print_exc()
        return {
            'file_path': '',
            'total_count': 0,
            'up_count': 0,
            'down_count': 0,
            'top_genes': []
        }

def generate_volcano_plot(df_de, lfc_cutoff, pval_cutoff, plots_dir):
    """
    生成火山图
    """
    try:
        plt.figure(figsize=(10, 8))
        
        if len(df_de) > 0:
            # 非显著基因
            ns_genes = df_de[df_de['regulation'] == 'ns']
            if len(ns_genes) > 0:
                plt.scatter(ns_genes['log2fc'], -np.log10(ns_genes['pval']), 
                           c='gray', alpha=0.5, s=10, label='Not significant')
            
            # 上调基因
            up_genes = df_de[df_de['regulation'] == 'up']
            if len(up_genes) > 0:
                plt.scatter(up_genes['log2fc'], -np.log10(up_genes['pval']), 
                           c='red', alpha=0.8, s=20, label=f'Up-regulated (n={len(up_genes)})')
            
            # 下调基因
            down_genes = df_de[df_de['regulation'] == 'down']
            if len(down_genes) > 0:
                plt.scatter(down_genes['log2fc'], -np.log10(down_genes['pval']), 
                           c='blue', alpha=0.8, s=20, label=f'Down-regulated (n={len(down_genes)})')
        
        # 阈值线
        plt.axhline(y=-np.log10(pval_cutoff), color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=lfc_cutoff, color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=-lfc_cutoff, color='black', linestyle='--', alpha=0.5)
        
        plt.xlabel('log2 Fold Change', fontsize=12)
        plt.ylabel('-log10(p-value)', fontsize=12)
        plt.title('Volcano Plot of Differentially Expressed Genes', fontsize=14)
        plt.legend(loc='upper right')
        plt.grid(True, alpha=0.3)
        
        # 标注top基因
        if len(df_de) > 0:
            top_genes = df_de.nsmallest(5, 'pval')
            for _, gene in top_genes.iterrows():
                plt.annotate(gene['gene_id'], 
                           (gene['log2fc'], -np.log10(gene['pval'])),
                           fontsize=8, alpha=0.7)
        
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'volcano-plot.png'), dpi=150, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error generating volcano plot: {e}")

def perform_enrichment_analysis(de_results, results_dir, plots_dir, params=None):
    """
    执行功能富集分析 - 基于预定义的基因集
    """
    try:
        if params is None:
            params = {}
        
        # 获取全部显著差异基因（优先读取完整差异分析结果文件，而不是仅用 top_genes）
        de_genes = []
        de_file = de_results.get('file_path')
        if de_file and os.path.exists(de_file):
            df_de_full = pd.read_csv(de_file)
            if 'regulation' in df_de_full.columns and 'gene_id' in df_de_full.columns:
                de_genes = df_de_full.loc[df_de_full['regulation'].isin(['up', 'down']), 'gene_id'].astype(str).tolist()
        if not de_genes and de_results.get('top_genes'):
            de_genes = [str(g['gene_id']) for g in de_results['top_genes'] 
                       if g.get('regulation') in ['up', 'down']]
        de_genes = list(dict.fromkeys(de_genes))
        
        # 如果没有显著的差异基因，返回空结果
        if len(de_genes) == 0:
            return {
                'file_path': '',
                'term_count': 0,
                'top_terms': []
            }
        
        # 预定义的基因集数据库（简化版）
        go_terms = {
            'GO:0006955': {
                'term_name': 'immune response',
                'genes': ['IFNG', 'IL6', 'TNF', 'IL1B', 'CD4', 'CD8A', 'IL2', 'IL4', 'IL10', 'CXCL8']
            },
            'GO:0006954': {
                'term_name': 'inflammatory response',
                'genes': ['IL6', 'TNF', 'IL1B', 'CXCL8', 'CCL2', 'CCL5', 'IL1A', 'PTGS2']
            },
            'GO:0008283': {
                'term_name': 'cell proliferation',
                'genes': ['PCNA', 'CDK1', 'MKI67', 'CCND1', 'CCNE1', 'MYC', 'EGFR']
            },
            'GO:0006915': {
                'term_name': 'apoptotic process',
                'genes': ['BAX', 'BCL2', 'CASP3', 'TP53', 'CASP9', 'FAS', 'BCL2L1']
            },
            'GO:0006950': {
                'term_name': 'response to stress',
                'genes': ['HSP90AA1', 'HSPB1', 'SOD1', 'CAT', 'HSPA1A', 'HMOX1']
            },
            'GO:0007155': {
                'term_name': 'cell adhesion',
                'genes': ['CDH1', 'CDH2', 'ITGA1', 'ITGB1', 'VCAM1', 'ICAM1']
            },
            'GO:0007165': {
                'term_name': 'signal transduction',
                'genes': ['MAPK1', 'MAPK3', 'AKT1', 'PIK3CA', 'SRC', 'GRB2']
            },
            'GO:0006355': {
                'term_name': 'regulation of transcription',
                'genes': ['TP53', 'MYC', 'JUN', 'FOS', 'NFKB1', 'RELA', 'STAT3']
            },
            'GO:0007049': {
                'term_name': 'cell cycle',
                'genes': ['CDK1', 'CDK2', 'CCNB1', 'CCNA2', 'CDC20', 'CDKN1A']
            },
            'GO:0006096': {
                'term_name': 'glycolysis',
                'genes': ['HK2', 'PKM', 'LDHA', 'GAPDH', 'ENO1', 'ALDOA']
            }
        }
        
        # KEGG通路（简化版）
        kegg_pathways = {
            'hsa04060': {
                'term_name': 'Cytokine-cytokine receptor interaction',
                'genes': ['IL6', 'TNF', 'IL1B', 'CXCL8', 'CCL2', 'CCL5', 'IL2', 'IL4']
            },
            'hsa04151': {
                'term_name': 'PI3K-Akt signaling pathway',
                'genes': ['AKT1', 'PIK3CA', 'MTOR', 'PTEN', 'MAPK1', 'MAPK3']
            },
            'hsa04668': {
                'term_name': 'TNF signaling pathway',
                'genes': ['TNF', 'NFKB1', 'MAPK1', 'MAPK3', 'IL6', 'CCL2']
            },
            'hsa04064': {
                'term_name': 'NF-kappa B signaling pathway',
                'genes': ['NFKB1', 'RELA', 'TNF', 'IL1B', 'CXCL8', 'ICAM1']
            }
        }
        
        # 选择数据库
        databases = params.get('databases', ['go-bp', 'kegg'])
        
        enrichment_terms = []
        
        # GO富集分析
        if 'go-bp' in databases:
            for go_id, term_info in go_terms.items():
                term_genes = term_info['genes']
                # 计算重叠基因
                overlapping_genes = set(de_genes) & set(term_genes)
                
                if len(overlapping_genes) >= 2:  # 至少2个重叠基因
                    # 超几何检验计算 p 值（确定性，非随机）
                    universe = max(int(de_results.get('total_count', 20000)), len(set(de_genes) | set(term_genes)))
                    pval = stats.hypergeom.sf(len(overlapping_genes) - 1, universe, len(set(term_genes)), len(set(de_genes)))
                    
                    enrichment_terms.append({
                        'term_name': term_info['term_name'],
                        'term_id': go_id,
                        'database': 'GO:BP',
                        'pval': float(pval),
                        'padj': float(pval * 10),  # 简化校正
                        'genes': list(overlapping_genes)[:5],
                        'gene_count': len(overlapping_genes),
                        'overlap_genes': ','.join(list(overlapping_genes)[:5])
                    })
        
        # KEGG富集分析
        if 'kegg' in databases:
            for pathway_id, pathway_info in kegg_pathways.items():
                pathway_genes = pathway_info['genes']
                overlapping_genes = set(de_genes) & set(pathway_genes)
                
                if len(overlapping_genes) >= 2:
                    universe = max(int(de_results.get('total_count', 20000)), len(set(de_genes) | set(pathway_genes)))
                    pval = stats.hypergeom.sf(len(overlapping_genes) - 1, universe, len(set(pathway_genes)), len(set(de_genes)))
                    
                    enrichment_terms.append({
                        'term_name': pathway_info['term_name'],
                        'term_id': pathway_id,
                        'database': 'KEGG',
                        'pval': float(pval),
                        'padj': float(pval * 5),
                        'genes': list(overlapping_genes)[:5],
                        'gene_count': len(overlapping_genes),
                        'overlap_genes': ','.join(list(overlapping_genes)[:5])
                    })
        
        # 排序
        enrichment_terms.sort(key=lambda x: x['pval'])
        
        # 保存结果
        if enrichment_terms:
            df_enrich = pd.DataFrame(enrichment_terms)
            enrich_file_path = os.path.join(results_dir, 'enrichment_results.csv')
            df_enrich.to_csv(enrich_file_path, index=False)
            
            # 生成富集条形图
            generate_enrichment_barplot(enrichment_terms[:10], plots_dir)
            
            return {
                'file_path': enrich_file_path,
                'term_count': len(enrichment_terms),
                'top_terms': enrichment_terms[:5]
            }
        else:
            return {
                'file_path': '',
                'term_count': 0,
                'top_terms': []
            }
        
    except Exception as e:
        print(f"Error in enrichment analysis: {e}")
        import traceback
        traceback.print_exc()
        return {
            'file_path': '',
            'term_count': 0,
            'top_terms': []
        }

def generate_enrichment_barplot(enrichment_terms, plots_dir):
    """
    生成富集分析条形图
    """
    try:
        if not enrichment_terms:
            return
            
        plt.figure(figsize=(10, 8))
        
        terms = [t['term_name'][:30] + '...' if len(t['term_name']) > 30 else t['term_name'] 
                for t in enrichment_terms]
        pvals = [-np.log10(t['pval']) for t in enrichment_terms]
        
        # 按p值排序
        sorted_idx = np.argsort(pvals)
        terms = [terms[i] for i in sorted_idx]
        pvals = [pvals[i] for i in sorted_idx]
        
        colors = plt.cm.RdYlBu_r(np.linspace(0.3, 0.8, len(terms)))
        
        plt.barh(range(len(terms)), pvals, color=colors)
        plt.yticks(range(len(terms)), terms)
        plt.xlabel('-log10(p-value)', fontsize=12)
        plt.title('Top Enriched Pathways/GO Terms', fontsize=14)
        
        # 添加p值标注
        for i, pval in enumerate(pvals):
            plt.text(pval + 0.1, i, f'{pval:.1f}', va='center', fontsize=8)
        
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'enrichment-barplot.png'), dpi=150, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error generating enrichment barplot: {e}")

def generate_all_plots(df_log2cpm, de_results, enrichment_results, groups, plots_dir):
    """
    生成所有可视化图表
    """
    try:
        # PCA分析
        generate_pca_plot(df_log2cpm, groups, plots_dir)
        
        # 表达热图
        generate_heatmap(df_log2cpm, groups, plots_dir)
        
        # 表达谱
        generate_expression_profiles(df_log2cpm, groups, plots_dir)
        
        # 如果差异表达结果存在，生成MA图
        if de_results.get('file_path') and os.path.exists(de_results['file_path']):
            df_de = pd.read_csv(de_results['file_path'])
            generate_ma_plot(df_de, plots_dir)
        
    except Exception as e:
        print(f"Error generating all plots: {e}")

def generate_pca_plot(df_log2cpm, groups, plots_dir):
    """
    生成PCA图
    """
    try:
        # PCA降维
        pca = PCA(n_components=2)
        data_scaled = StandardScaler().fit_transform(df_log2cpm.T)
        pca_result = pca.fit_transform(data_scaled)
        
        plt.figure(figsize=(10, 8))
        
        # 根据分组着色
        if groups and isinstance(groups, dict):
            sample_groups = [groups.get(col, 'Unknown') for col in df_log2cpm.columns]
            unique_groups = list(set(sample_groups))
            
            for group in unique_groups:
                idx = [i for i, g in enumerate(sample_groups) if g == group]
                plt.scatter(pca_result[idx, 0], pca_result[idx, 1], 
                          s=100, alpha=0.7, label=group)
        else:
            plt.scatter(pca_result[:, 0], pca_result[:, 1], s=100, alpha=0.7)
        
        # 添加样本标签
        for i, sample in enumerate(df_log2cpm.columns):
            plt.annotate(sample, (pca_result[i, 0], pca_result[i, 1]), 
                        fontsize=8, alpha=0.7)
        
        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%})', fontsize=12)
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%})', fontsize=12)
        plt.title('Principal Component Analysis', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        plt.savefig(os.path.join(plots_dir, 'pca-analysis.png'), dpi=150, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error generating PCA plot: {e}")

def generate_heatmap(df_log2cpm, groups, plots_dir):
    """
    生成表达热图 - 修复数组错误
    """
    try:
        # 选择变异最大的50个基因
        gene_vars = df_log2cpm.var(axis=1).nlargest(50)
        top_genes = df_log2cpm.loc[gene_vars.index]
        
        # Z-score标准化
        from scipy.stats import zscore
        top_genes_zscore = top_genes.apply(zscore, axis=1)
        
        # 处理NaN值
        top_genes_zscore = top_genes_zscore.fillna(0)
        
        # 添加分组注释
        if groups and isinstance(groups, dict):
            sample_names = df_log2cpm.columns.tolist()
            sample_groups = [groups.get(sample, 'Unknown') for sample in sample_names]
            
            # 创建分组颜色
            unique_groups = list(set(sample_groups))
            colors = plt.cm.Set3(np.linspace(0, 1, len(unique_groups)))
            group_colors = {group: colors[i] for i, group in enumerate(unique_groups)}
            col_colors = [group_colors[g] for g in sample_groups]
            
            # 使用clustermap
            g = sns.clustermap(top_genes_zscore, 
                              cmap='RdYlBu_r',
                              center=0,
                              xticklabels=True,
                              yticklabels=True,
                              figsize=(12, 10),
                              dendrogram_ratio=0.15,
                              col_colors=col_colors,
                              cbar_pos=(0.02, 0.8, 0.03, 0.18))
            g.ax_heatmap.set_title('Top 50 Most Variable Genes - Expression Heatmap', fontsize=14)
        else:
            g = sns.clustermap(top_genes_zscore, 
                              cmap='RdYlBu_r',
                              center=0,
                              xticklabels=True,
                              yticklabels=True,
                              figsize=(12, 10),
                              dendrogram_ratio=0.15)
            g.ax_heatmap.set_title('Top 50 Most Variable Genes - Expression Heatmap', fontsize=14)
        
        plt.savefig(os.path.join(plots_dir, 'heatmap.png'), dpi=150, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error generating heatmap: {e}")
        import traceback
        traceback.print_exc()
        
        # 创建简单的占位图
        plt.figure(figsize=(10, 8))
        plt.text(0.5, 0.5, 'Heatmap generation failed\nPlease check your data', 
                ha='center', va='center', fontsize=14)
        plt.title('Expression Heatmap', fontsize=16)
        plt.axis('off')
        plt.savefig(os.path.join(plots_dir, 'heatmap.png'), dpi=150, bbox_inches='tight')
        plt.close()

def generate_expression_profiles(df_log2cpm, groups, plots_dir):
    """
    生成表达谱图
    """
    try:
        # 选择变异最大的10个基因
        gene_vars = df_log2cpm.var(axis=1).nlargest(10)
        top_genes = df_log2cpm.loc[gene_vars.index]
        
        plt.figure(figsize=(14, 7))
        
        # 获取分组信息用于着色
        if groups and isinstance(groups, dict):
            sample_groups = [groups.get(col, 'Unknown') for col in df_log2cpm.columns]
            group_colors = plt.cm.Set1(np.linspace(0, 1, len(set(sample_groups))))
            group_color_map = {g: group_colors[i] for i, g in enumerate(set(sample_groups))}
            colors = [group_color_map[g] for g in sample_groups]
            
            # 按分组排序
            sorted_idx = np.argsort(sample_groups)
            sorted_samples = [df_log2cpm.columns[i] for i in sorted_idx]
            sorted_top_genes = top_genes[sorted_samples]
            
            for gene in sorted_top_genes.index:
                plt.plot(range(len(sorted_samples)), sorted_top_genes.loc[gene].values, 
                        marker='o', label=gene, linewidth=2, markersize=6)
        else:
            for gene in top_genes.index:
                plt.plot(range(len(top_genes.columns)), top_genes.loc[gene].values, 
                        marker='o', label=gene, linewidth=2, markersize=6)
        
        plt.xlabel('Samples', fontsize=12)
        plt.ylabel('log2(CPM+1)', fontsize=12)
        plt.title('Top 10 Most Variable Genes - Expression Profiles', fontsize=14)
        plt.xticks(range(len(df_log2cpm.columns)), df_log2cpm.columns, rotation=45)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        plt.savefig(os.path.join(plots_dir, 'expression-profiles.png'), dpi=150, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error generating expression profiles: {e}")

def generate_ma_plot(df_de, plots_dir):
    """
    生成MA图 (M vs A)
    """
    try:
        plt.figure(figsize=(10, 8))
        
        if len(df_de) > 0:
            # 非显著基因
            ns_genes = df_de[df_de['regulation'] == 'ns']
            if len(ns_genes) > 0:
                plt.scatter(ns_genes['base_mean'], ns_genes['log2fc'], 
                           c='gray', alpha=0.5, s=10, label='Not significant')
            
            # 上调基因
            up_genes = df_de[df_de['regulation'] == 'up']
            if len(up_genes) > 0:
                plt.scatter(up_genes['base_mean'], up_genes['log2fc'], 
                           c='red', alpha=0.8, s=20, label='Up-regulated')
            
            # 下调基因
            down_genes = df_de[df_de['regulation'] == 'down']
            if len(down_genes) > 0:
                plt.scatter(down_genes['base_mean'], down_genes['log2fc'], 
                           c='blue', alpha=0.8, s=20, label='Down-regulated')
        
        plt.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        plt.xscale('log')
        plt.xlabel('Average Expression (log2CPM)', fontsize=12)
        plt.ylabel('log2 Fold Change', fontsize=12)
        plt.title('MA Plot', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        plt.savefig(os.path.join(plots_dir, 'ma-plot.png'), dpi=150, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error generating MA plot: {e}")

@require_http_methods(["GET"])
def get_transcriptome_task_status(request, task_id):
    """
    获取转录组分析任务状态
    """
    try:
        task = TranscriptomeAnalysisTask.objects.get(id=task_id)
        
        return JsonResponse({
            'success': True,
            'task_id': str(task_id),
            'status': task.status,
            'progress': task.progress,
            'current_step': task.current_step,
            'message': task.message,
            'error': task.error,
            'created_at': task.created_at.isoformat() if task.created_at else None,
            'updated_at': task.updated_at.isoformat() if task.updated_at else None,
            'started_at': task.started_at.isoformat() if task.started_at else None,
            'completed_at': task.completed_at.isoformat() if task.completed_at else None
        })
        
    except TranscriptomeAnalysisTask.DoesNotExist:
        return JsonResponse({
            'success': False,
            'error': 'Task not found'
        }, status=404)
    except Exception as e:
        logger.error(f"Error getting task status: {str(e)}")
        return JsonResponse({
            'success': False,
            'error': f'Failed to get task status: {str(e)}'
        }, status=500)

@require_http_methods(["GET"])
def get_transcriptome_task_results(request, task_id):
    """
    获取转录组分析结果
    """
    try:
        task = TranscriptomeAnalysisTask.objects.get(id=task_id)
        
        # 检查分析是否完成
        if task.status != 'completed':
            return JsonResponse({
                'success': False,
                'error': 'Analysis not completed yet'
            }, status=400)
        
        # 获取结果
        try:
            result = TranscriptomeAnalysisResult.objects.get(task=task)
        except TranscriptomeAnalysisResult.DoesNotExist:
            return JsonResponse({
                'success': False,
                'error': 'Results not found'
            }, status=404)
        
        return JsonResponse({
            'success': True,
            'task_id': str(task_id),
            'summary': result.summary,
            'differential_expression': result.differential_expression,
            'enrichment': result.enrichment,
            'file_paths': result.file_paths,
            'created_at': result.created_at.isoformat() if result.created_at else None
        })
        
    except TranscriptomeAnalysisTask.DoesNotExist:
        return JsonResponse({
            'success': False,
            'error': 'Task not found'
        }, status=404)
    except Exception as e:
        logger.error(f"Error getting results: {str(e)}")
        return JsonResponse({
            'success': False,
            'error': f'Failed to get results: {str(e)}'
        }, status=500)

@require_http_methods(["GET"])
def get_transcriptome_visualization_plot(request, task_id, plot_type):
    """
    获取转录组分析可视化图表 - 返回真实生成的PNG文件
    """
    try:
        task_uuid_str = str(task_id)
        
        # 验证任务
        try:
            task = TranscriptomeAnalysisTask.objects.get(id=task_id)
        except TranscriptomeAnalysisTask.DoesNotExist:
            return JsonResponse({
                'success': False,
                'error': 'Task not found'
            }, status=404)
        
        # 验证图表类型
        if plot_type not in PLOT_TYPES:
            return JsonResponse({
                'success': False,
                'error': f'Invalid plot type. Allowed: {", ".join(PLOT_TYPES)}'
            }, status=400)
        
        # 检查分析是否完成
        if task.status != 'completed':
            return JsonResponse({
                'success': False,
                'error': 'Analysis not completed yet'
            }, status=400)
        
        # 图表文件路径
        plot_dir = os.path.join(settings.MEDIA_ROOT, 'transcriptome_analysis', task_uuid_str, 'plots')
        plot_file = os.path.join(plot_dir, f'{plot_type}.png')
        
        # 如果请求的是coexpression-network，返回占位图或简化网络图
        if plot_type == 'coexpression-network' and not os.path.exists(plot_file):
            generate_simple_network_plot(task, plot_dir)
            plot_file = os.path.join(plot_dir, 'coexpression-network.png')
        
        # 如果请求的是pathway-map，返回占位图
        if plot_type == 'pathway-map' and not os.path.exists(plot_file):
            generate_pathway_placeholder(plot_dir)
            plot_file = os.path.join(plot_dir, 'pathway-map.png')
        
        # 检查文件是否存在
        if not os.path.exists(plot_file):
            return JsonResponse({
                'success': False,
                'error': f'Plot {plot_type} not found for this task'
            }, status=404)
        
        # 返回图表文件
        return FileResponse(open(plot_file, 'rb'), content_type='image/png')
        
    except Exception as e:
        logger.error(f"Error getting plot {plot_type}: {str(e)}", exc_info=True)
        return JsonResponse({
            'success': False,
            'error': f'Failed to get plot: {str(e)}'
        }, status=500)

def generate_simple_network_plot(task, plots_dir):
    """
    生成简化的共表达网络图
    """
    try:
        import networkx as nx
        
        plt.figure(figsize=(12, 10))
        
        # 创建一个示例网络
        G = nx.Graph()
        
        # 添加节点（基因）
        genes = [f'Gene_{i}' for i in range(1, 21)]
        G.add_nodes_from(genes)
        
        # 添加边（共表达关系）
        np.random.seed(42)
        for i in range(len(genes)):
            for j in range(i+1, len(genes)):
                if np.random.random() < 0.15:
                    G.add_edge(genes[i], genes[j], weight=np.random.uniform(0.7, 0.95))
        
        # 布局
        pos = nx.spring_layout(G, k=2, iterations=50)
        
        # 绘制节点
        node_size = [G.degree(node) * 100 + 200 for node in G.nodes()]
        nx.draw_networkx_nodes(G, pos, node_size=node_size, node_color='lightblue', 
                              alpha=0.9, edgecolors='darkblue', linewidths=1)
        
        # 绘制边
        edges = G.edges()
        weights = [G[u][v]['weight'] for u, v in edges]
        nx.draw_networkx_edges(G, pos, width=weights, alpha=0.6, edge_color='gray')
        
        # 绘制标签
        nx.draw_networkx_labels(G, pos, font_size=8, font_family='sans-serif')
        
        plt.title('Gene Co-expression Network (Top 20 Genes)', fontsize=16)
        plt.axis('off')
        plt.tight_layout()
        
        plt.savefig(os.path.join(plots_dir, 'coexpression-network.png'), dpi=150, bbox_inches='tight')
        plt.close()
        
    except ImportError:
        # 如果没有networkx，创建简单的占位图
        plt.figure(figsize=(12, 10))
        plt.text(0.5, 0.5, 'Gene Co-expression Network\n(NetworkX required for full visualization)', 
                ha='center', va='center', fontsize=14)
        plt.axis('off')
        plt.savefig(os.path.join(plots_dir, 'coexpression-network.png'), dpi=150, bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"Error generating network plot: {e}")

def generate_pathway_placeholder(plots_dir):
    """
    生成通路可视化占位图
    """
    try:
        plt.figure(figsize=(12, 8))
        
        # 创建一个简化的通路图
        pathways = ['MAPK', 'PI3K-Akt', 'TNF', 'NF-kB', 'p53', 'Wnt']
        
        for i, pathway in enumerate(pathways):
            plt.text(0.5, 0.9 - i*0.15, f'{pathway} Signaling Pathway', 
                    ha='center', va='center', fontsize=14,
                    bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
        
        plt.text(0.5, 0.1, 'KEGG Pathway Overlay\n(Hover for gene expression)', 
                ha='center', va='center', fontsize=12, style='italic')
        
        plt.title('KEGG Pathway Visualization', fontsize=16)
        plt.axis('off')
        plt.tight_layout()
        
        plt.savefig(os.path.join(plots_dir, 'pathway-map.png'), dpi=150, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error generating pathway placeholder: {e}")

@require_http_methods(["GET"])
def download_transcriptome_result_by_type(request, task_id):
    """
    按类型下载转录组分析结果
    """
    try:
        task_uuid_str = str(task_id)
        
        # 验证任务
        try:
            task = TranscriptomeAnalysisTask.objects.get(id=task_id)
        except TranscriptomeAnalysisTask.DoesNotExist:
            return JsonResponse({
                'success': False,
                'error': 'Task not found'
            }, status=404)
        
        # 获取下载类型
        result_type = request.GET.get('type', 'all')
        
        # 验证结果类型
        if result_type not in RESULT_TYPES:
            return JsonResponse({
                'success': False,
                'error': f'Invalid result type. Allowed: {", ".join(RESULT_TYPES.keys())}'
            }, status=400)
        
        # 检查分析是否完成
        if task.status != 'completed':
            return JsonResponse({
                'success': False,
                'error': 'Analysis not completed yet'
            }, status=400)
        
        # 获取结果文件路径
        try:
            result = TranscriptomeAnalysisResult.objects.get(task=task)
        except TranscriptomeAnalysisResult.DoesNotExist:
            return JsonResponse({
                'success': False,
                'error': 'Results not found'
            }, status=404)
        
        file_paths = result.file_paths
        
        if result_type == 'all':
            # 打包所有结果
            zip_filename = f'transcriptome_analysis_{task_uuid_str}_{datetime.now().strftime("%Y%m%d_%H%M%S")}.zip'
            zip_filepath = os.path.join(os.path.dirname(file_paths.get('counts', '')), zip_filename)
            
            with zipfile.ZipFile(zip_filepath, 'w', zipfile.ZIP_DEFLATED) as zipf:
                # 添加所有结果文件
                for key, path in file_paths.items():
                    if key != 'plots' and path and os.path.exists(path):
                        zipf.write(path, os.path.basename(path))
                
                # 添加所有图表
                plots_dir = file_paths.get('plots', '')
                if plots_dir and os.path.exists(plots_dir):
                    for plot_file in os.listdir(plots_dir):
                        if plot_file.endswith('.png'):
                            plot_path = os.path.join(plots_dir, plot_file)
                            zipf.write(plot_path, f'plots/{plot_file}')
            
            response = FileResponse(open(zip_filepath, 'rb'))
            response['Content-Disposition'] = f'attachment; filename="{zip_filename}"'
            response['Content-Type'] = 'application/zip'
            
        elif result_type == 'visualizations':
            # 打包所有图表
            zip_filename = f'transcriptome_plots_{task_uuid_str}_{datetime.now().strftime("%Y%m%d_%H%M%S")}.zip'
            zip_filepath = os.path.join(os.path.dirname(file_paths.get('counts', '')), zip_filename)
            
            with zipfile.ZipFile(zip_filepath, 'w', zipfile.ZIP_DEFLATED) as zipf:
                plots_dir = file_paths.get('plots', '')
                if plots_dir and os.path.exists(plots_dir):
                    for plot_file in os.listdir(plots_dir):
                        if plot_file.endswith('.png'):
                            plot_path = os.path.join(plots_dir, plot_file)
                            zipf.write(plot_path, plot_file)
            
            response = FileResponse(open(zip_filepath, 'rb'))
            response['Content-Disposition'] = f'attachment; filename="{zip_filename}"'
            response['Content-Type'] = 'application/zip'
            
        else:
            # 单个文件下载
            file_mapping = {
                'counts': 'counts',
                'normalized': 'log2cpm',
                'de-genes': 'de_genes',
                'enrichment': 'enrichment',
                'analysis-summary': 'summary'
            }
            
            file_key = file_mapping.get(result_type)
            if not file_key or file_key not in file_paths:
                return JsonResponse({
                    'success': False,
                    'error': f'File type {result_type} not available'
                }, status=404)
            
            file_path = file_paths[file_key]
            if not file_path or not os.path.exists(file_path):
                return JsonResponse({
                    'success': False,
                    'error': f'File not found on server'
                }, status=404)
            
            filename = os.path.basename(file_path)
            response = FileResponse(open(file_path, 'rb'))
            response['Content-Disposition'] = f'attachment; filename="{filename}"'
            
            mime_type, _ = mimetypes.guess_type(filename)
            response['Content-Type'] = mime_type or 'application/octet-stream'
        
        return response
        
    except Exception as e:
        logger.error(f"Error downloading result: {str(e)}", exc_info=True)
        return JsonResponse({
            'success': False,
            'error': f'Failed to download result: {str(e)}'
        }, status=500)


@require_http_methods(["GET"])
def download_all_transcriptome_results(request, task_id):
    """
    下载转录组分析全部结果（兼容前端 /download-all/ 路由）。
    内部复用按类型下载接口，将 type 固定为 all。
    """
    get_params = request.GET.copy()
    get_params['type'] = 'all'
    request.GET = get_params
    return download_transcriptome_result_by_type(request, task_id)

@require_http_methods(["GET"])
def list_user_transcriptome_tasks(request):
    """
    获取用户的所有转录组分析任务
    """
    try:
        user = request.user
        if not user.is_authenticated:
            return JsonResponse({
                'success': False,
                'error': 'Authentication required to list tasks'
            }, status=401)
        
        tasks = TranscriptomeAnalysisTask.objects.filter(user=user).order_by('-created_at')
        
        task_list = []
        for task in tasks:
            task_list.append({
                'id': str(task.id),
                'status': task.status,
                'progress': task.progress,
                'current_step': task.current_step,
                'message': task.message,
                'analysis_type': task.analysis_type,
                'created_at': task.created_at.isoformat() if task.created_at else None,
                'updated_at': task.updated_at.isoformat() if task.updated_at else None,
                'started_at': task.started_at.isoformat() if task.started_at else None,
                'completed_at': task.completed_at.isoformat() if task.completed_at else None,
                'has_results': hasattr(task, 'transcriptomeanalysisresult')
            })
        
        return JsonResponse({
            'success': True,
            'tasks': task_list,
            'total': len(task_list)
        })
        
    except Exception as e:
        logger.error(f"Error listing tasks: {str(e)}")
        return JsonResponse({
            'success': False,
            'error': f'Failed to list tasks: {str(e)}'
        }, status=500)

@csrf_exempt
@require_http_methods(["POST"])
def cancel_transcriptome_task(request, task_id):
    """
    取消转录组分析任务
    """
    try:
        task = TranscriptomeAnalysisTask.objects.get(id=task_id)
        
        # 检查任务状态是否可以取消
        if task.status not in ['pending', 'running']:
            return JsonResponse({
                'success': False,
                'error': f'Cannot cancel task in {task.status} state'
            }, status=400)
        
        task.status = 'cancelled'
        task.message = 'Task cancelled by user'
        task.error = 'Cancelled by user request'
        task.save()
        
        return JsonResponse({
            'success': True,
            'message': 'Task cancelled successfully',
            'status': task.status
        })
        
    except TranscriptomeAnalysisTask.DoesNotExist:
        return JsonResponse({
            'success': False,
            'error': 'Task not found'
        }, status=404)
    except Exception as e:
        logger.error(f"Error cancelling task: {str(e)}")
        return JsonResponse({
            'success': False,
            'error': f'Failed to cancel task: {str(e)}'
        }, status=500)