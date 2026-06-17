import os
import subprocess
import uuid
import json
import logging
from datetime import datetime # 可能不需要，用 Django 的 timezone
from django.shortcuts import render
from django.http import JsonResponse, FileResponse, Http404
from django.conf import settings
from django.core.files.storage import default_storage
from django.views.decorators.csrf import csrf_exempt
from django.utils.encoding import escape_uri_path
from django.utils import timezone # Import Django's timezone utility
from io import BytesIO
import matplotlib.pyplot as plt
import pandas as pd # 假设您已安装pandas用于处理CSV结果
import seaborn as sns # 假设您已安装seaborn用于绘图
from Bio import SeqIO # 引入Biopython用于处理FASTA文件
import zipfile # 用于打包下载

# --- 导入 BlastTask 模型 ---
from ..models import BlastTask # 假设模型在同级的 models.py 文件中

logger = logging.getLogger(__name__)

# --- 配置 ---
# 请根据您的实际环境修改这些路径
BLAST_DATABASE_ROOT = "/www/wwwroot/default/animal_hormone/hormone_app/blast"
BLAST_PROGRAMS_PATH = "/www/wwwroot/default/animal_hormone/hormone_app/blast/bin" 
RESULTS_STORAGE_PATH = os.path.join(settings.MEDIA_ROOT, 'blast_results')
# 自动创建目录（确保存在）
os.makedirs(RESULTS_STORAGE_PATH, exist_ok=True)

# --- 数据库映射 ---
# 将前端发送的数据库名称映射到实际的数据库文件夹名称
# 现在只需要映射到数据库目录名，因为我们会用 blastdbcmd 从数据库本身提取序列

# --- 数据库映射 - 使用完整的数据库前缀路径 ---
DATABASE_MAPPING = {
    # 核酸数据库 - 使用别名文件（已修复）
    'nt': os.path.join(BLAST_DATABASE_ROOT, 'animal_hormone_db', 'GCF_000001405.39_top_level'),
    'refseq_rna': os.path.join(BLAST_DATABASE_ROOT, 'animal_hormone_db', 'GCF_000001405.39_top_level'),
    'custom': os.path.join(BLAST_DATABASE_ROOT, 'animal_hormone_db', 'GCF_000001405.39_top_level'),
    
    # 蛋白质数据库
    'nr': os.path.join(BLAST_DATABASE_ROOT, 'animal_hormone_p', 'human_hormones_blastdb'),
    'refseq_protein': os.path.join(BLAST_DATABASE_ROOT, 'animal_hormone_p', 'human_full_sprot_blastdb'),
    'swissprot': os.path.join(BLAST_DATABASE_ROOT, 'animal_hormone_p', 'human_hormones_blastdb'),
}

# 数据库类型验证常量
NUCLEOTIDE_DATABASES = {'nt', 'refseq_rna', 'custom'}
PROTEIN_DATABASES = {'nr', 'refseq_protein', 'swissprot'}

# 程序类型验证常量
NUCLEOTIDE_PROGRAMS = {'blastn'}
PROTEIN_PROGRAMS = {'blastp', 'blastx', 'tblastn', 'tblastx'}

# --- 移除模拟数据库存储 ---
# TASKS_DB = {} # 删除这一行

def run_blast_command(query_file_path, db_name, output_file_path, task_params, outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'):
    """生产级 BLAST 命令执行：校验程序/数据库/参数，失败返回明确日志。"""
    try:
        db_prefix_path = DATABASE_MAPPING.get(db_name)
        if not db_prefix_path:
            return False, f"Database mapping not found: {db_name}"
        program = str(task_params.get('program', 'blastn'))
        program_path = os.path.join(BLAST_PROGRAMS_PATH, program)
        if program not in NUCLEOTIDE_PROGRAMS.union(PROTEIN_PROGRAMS):
            return False, f"Unsupported BLAST program: {program}"
        if not os.path.exists(program_path) or not os.access(program_path, os.X_OK):
            return False, f"BLAST executable not found or not executable: {program_path}"
        try:
            evalue_float = float(task_params.get('evalue', '10.0'))
            if evalue_float <= 0:
                raise ValueError
        except Exception:
            return False, "Invalid evalue; must be a positive number"
        try:
            max_hits_int = int(task_params.get('max_hits', '50'))
            if max_hits_int < 1 or max_hits_int > 5000:
                raise ValueError
        except Exception:
            return False, "Invalid max_hits; must be integer 1-5000"
        cmd = [program_path, '-query', query_file_path, '-db', db_prefix_path, '-out', output_file_path,
               '-outfmt', str(outfmt), '-evalue', str(evalue_float), '-max_target_seqs', str(max_hits_int)]
        logger.info(f"Running BLAST command: {' '.join(cmd)}")
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=1800)
        if result.stderr:
            logger.info(f"BLAST stderr: {result.stderr[:2000]}")
        return True, ''
    except subprocess.CalledProcessError as e:
        msg = e.stderr or e.stdout or str(e)
        logger.error(f"BLAST command failed: {msg}")
        return False, msg
    except Exception as e:
        logger.error(f"An error occurred while running BLAST: {str(e)}")
        return False, str(e)

def fetch_fasta_hits_via_blastdbcmd(hit_ids, db_path):
    """使用 blastdbcmd 从数据库中提取对应的FASTA序列"""
    hits = []
    if not hit_ids:
        logger.warning("No hit IDs provided to fetch_fasta_hits_via_blastdbcmd.")
        return hits

    try:
        # 使用 -target_only 确保只提取指定的序列
        # 使用默认的 outfmt 输出完整FASTA格式
        cmd = [
            os.path.join(BLAST_PROGRAMS_PATH, 'blastdbcmd'),
            '-db', db_path,
            '-entry_batch', '-',
            '-target_only'  # 只提取指定的序列
        ]
        
        id_input = '\n'.join(hit_ids) + '\n'
        logger.info(f"Running blastdbcmd: {' '.join(cmd)}")

        result = subprocess.run(cmd, input=id_input, text=True, capture_output=True, check=True)
        
        if not result.stdout.strip():
            logger.warning("blastdbcmd returned empty output.")
            return []
            
        # 使用 StringIO 解析FASTA内容
        from io import StringIO
        hits = list(SeqIO.parse(StringIO(result.stdout), "fasta"))
        
        logger.info(f"Fetched {len(hits)} sequences using blastdbcmd.")
        return hits

    except subprocess.CalledProcessError as e:
        logger.error(f"blastdbcmd failed for database {db_path}: {e.stderr}")
        # 尝试更详细的错误信息
        logger.error(f"Command output: {e.stdout}")
        return []
    except Exception as e:
        logger.error(f"Error fetching FASTA hits via blastdbcmd for database {db_path}: {str(e)}")
        return []


@csrf_exempt
def blast_create_task(request):
    """
    创建一个新的BLAST分析任务。
    此函数现在只接收任务的基本信息（名称、描述、序列类型），用于初始化一个任务。
    序列和数据库等详细参数将在 blast_run_analysis 时提供。
    """
    logger.info(f"[CREATE_TASK] Received POST request to create task. Method: {request.method}")

    if request.method != 'POST':
        logger.warning(f"[CREATE_TASK] Invalid method {request.method} for create task endpoint.")
        return JsonResponse({'error': 'Method not allowed'}, status=405)

    try:
        # --- 修改：使用 request.POST.get() 替代 json.loads(request.body) ---
        # 因为前端发送的是 application/x-www-form-urlencoded 格式
        task_name = request.POST.get('task_name', f'BLAST_Task_{uuid.uuid4().hex[:8]}')
        description = request.POST.get('description', '')
        sequence_type = request.POST.get('sequence_type', 'protein') # 'protein' 或 'nucleotide'

        # 验证输入
        if not task_name.strip():
            logger.error("[CREATE_TASK] Task name cannot be empty.")
            return JsonResponse({'error': 'Task name cannot be empty.'}, status=400)

        task_id = str(uuid.uuid4())
        
        # --- 修改：不再使用内存字典，而是创建数据库记录 ---
        task_obj = BlastTask.objects.create(
            id=task_id,
            name=task_name,
            description=description,
            sequence_type=sequence_type,
            # status 默认为 'created'
            # 其他字段留空
        )
        # --- 修改结束 ---

        logger.info(f"[CREATE_TASK] Created new BLAST task object in DB: {task_id}, name: {task_name}, type: {sequence_type}")
        # 不再记录内存字典状态
        
        # 返回包含 task_id 的 JSON 响应
        return JsonResponse({'task_id': task_id, 'message': 'Task created successfully'}, status=201)

    # --- 修改：移除 json.JSONDecodeError 处理，因为不再需要 ---
    except Exception as e:
        logger.error(f"[CREATE_TASK] Error creating task: {str(e)}")
        return JsonResponse({'error': 'Internal server error'}, status=500)


@csrf_exempt
def blast_run_analysis(request, task_id):
    """
    为指定任务提交序列和参数，并执行 BLAST 分析。
    这个函数现在接收序列、数据库等参数，并执行分析。
    """
    logger.info(f"[RUN_ANALYSIS] Received POST request to run analysis for task ID: {task_id}. Method: {request.method}")

    if request.method != 'POST':
        logger.warning(f"[RUN_ANALYSIS] Invalid method {request.method} for run analysis endpoint.")
        return JsonResponse({'error': 'Method not allowed'}, status=405)

    # --- 修改：从数据库查询任务对象 ---
    try:
        task_obj = BlastTask.objects.get(id=task_id)
    except BlastTask.DoesNotExist:
        logger.error(f"[RUN_ANALYSIS] Task ID '{task_id}' not found in BlastTask model.")
        return JsonResponse({'error': 'Task not found'}, status=404)
    # --- 修改结束 ---

    if task_obj.status in ['running', 'completed']:
        logger.warning(f"[RUN_ANALYSIS] Attempted to run task '{task_id}' which is already in '{task_obj.status}' state.")
        return JsonResponse({'error': 'Task is already running or completed'}, status=400)

    try:
        # --- 在 run_analysis 阶段接收参数 ---
        # 从前端请求体获取完整参数 (假设前端此时发送的是JSON)
        data = json.loads(request.body)
        sequence = data.get('sequence', '')
        requested_database = data.get('database', 'nt')
        program = data.get('program', 'blastn')
        evalue = data.get('evalue', '10.0')
        max_hits = data.get('max_hits', '50')

        logger.info(f"[RUN_ANALYSIS] Received parameters for task '{task_id}': database={requested_database}, program={program}, evalue={evalue}, max_hits={max_hits}, sequence_length={len(sequence)}")
        logger.info(f"[RUN_ANALYSIS] BLAST program: {program}, E-value: {evalue}, Max hits: {max_hits}")

        # --- 新增：检测序列内容类型并记录警告 ---
        def _detect_sequence_type(seq_str):
            if not seq_str:
                return None
            seq_upper = seq_str.upper().replace('\n', '').replace('\r', '').replace(' ', '')
            if not seq_upper:
                return None
            protein_chars = set('EFILPQXZ*')
            potential_nucleo_chars = set('ATGCUNRYKMSWBDHVatgcunrykmswbdhv')
            char_set = set(seq_upper)
            if char_set.intersection(protein_chars):
                return 'protein'
            non_nucleo = char_set - potential_nucleo_chars
            nucleo_count = sum(seq_upper.count(c) for c in 'ATGCUatgcu')
            total_len = len(seq_upper)
            if total_len > 0:
                ratio = nucleo_count / total_len
                if ratio > 0.85:
                    return 'nucleotide'
            return 'protein'

        detected_seq_type = _detect_sequence_type(sequence)
        if detected_seq_type and detected_seq_type != task_obj.sequence_type:
             logger.warning(f"[RUN_ANALYSIS] For task '{task_id}', user-specified sequence type was '{task_obj.sequence_type}', "
                            f"but auto-detected sequence content suggests type '{detected_seq_type}'. "
                            f"The submitted sequence is '{sequence[:100]}...'. "
                            f"This might lead to unexpected analysis results.")
        # --- END 新增 ---

        # 验证数据库名称是否在映射中
        if requested_database not in DATABASE_MAPPING:
             logger.error(f"[RUN_ANALYSIS] Unsupported database requested: {requested_database}")
             return JsonResponse({'error': f'Unsupported database: {requested_database}'}, status=400)

        # 获取实际使用的数据库路径（前缀路径）
        actual_database_path = DATABASE_MAPPING[requested_database]
        logger.info(f"[RUN_ANALYSIS] Actual database path: {actual_database_path}")

        # --- 验证 program 与 task 的 sequence_type 的兼容性 ---
        if task_obj.sequence_type == 'nucleotide':
            if program not in NUCLEOTIDE_PROGRAMS:
                logger.error(f"[RUN_ANALYSIS] Incompatible program '{program}' for nucleotide task type on task '{task_id}'.")
                return JsonResponse({'error': f'Protein search program ({program}) cannot be used with nucleotide sequence type.'}, status=400)
        else: # protein sequence type
            if program not in PROTEIN_PROGRAMS:
                logger.error(f"[RUN_ANALYSIS] Incompatible program '{program}' for protein task type on task '{task_id}'.")
                return JsonResponse({'error': f'Nucleotide search program ({program}) cannot be used with protein sequence type.'}, status=400)

        # --- 验证 program 与 database 类型的兼容性 ---
        if requested_database in NUCLEOTIDE_DATABASES:
            if program not in NUCLEOTIDE_PROGRAMS:
                logger.error(f"[RUN_ANALYSIS] Incompatible program '{program}' for nucleotide database '{requested_database}' on task '{task_id}'.")
                return JsonResponse({'error': f'Protein search program ({program}) cannot be used with nucleotide database ({requested_database}).'}, status=400)
        elif requested_database in PROTEIN_DATABASES:
            if program not in PROTEIN_PROGRAMS:
                logger.error(f"[RUN_ANALYSIS] Incompatible program '{program}' for protein database '{requested_database}' on task '{task_id}'.")
                return JsonResponse({'error': f'Nucleotide search program ({program}) cannot be used with protein database ({requested_database}).'}, status=400)
        else:
            logger.error(f"[RUN_ANALYSIS] Unknown database type for '{requested_database}'")
            return JsonResponse({'error': f'Unknown database type: {requested_database}'}, status=400)

        if not sequence:
            logger.error(f"[RUN_ANALYSIS] Sequence is required for task '{task_id}'.")
            return JsonResponse({'error': 'Sequence is required'}, status=400)

        # 更新任务状态
        task_obj.status = 'running'
        
        # 保存序列到临时文件
        query_filename = f"{task_id}_query.fasta"
        query_file_path = os.path.join(RESULTS_STORAGE_PATH, query_filename)
        with open(query_file_path, 'w') as f:
            f.write(f">query_seq\n{sequence}\n")

        # 设置任务参数到数据库对象
        task_obj.sequence_file = query_file_path
        task_obj.database_requested = requested_database
        task_obj.database_actual = actual_database_path  # 保存完整的数据库前缀路径
        task_obj.program = program
        task_obj.evalue = evalue
        task_obj.max_hits = max_hits

        logger.info(f"[RUN_ANALYSIS] Updated task '{task_id}' status to 'running' and saved parameters to DB.")

        # 准备输出文件路径
        result_file_path = os.path.join(RESULTS_STORAGE_PATH, f"{task_id}_raw.out")
        csv_result_path = os.path.join(RESULTS_STORAGE_PATH, f"{task_id}_results.csv")
        xml_report_path = os.path.join(RESULTS_STORAGE_PATH, f"{task_id}_report.xml")
        fasta_hits_path = os.path.join(RESULTS_STORAGE_PATH, f"{task_id}_hits.fasta")

        # 运行 -outfmt 6 获取TSV结果
        success_tsv, blast_error = run_blast_command(
            task_obj.sequence_file, 
            requested_database,
            result_file_path,
            {'program': task_obj.program, 'evalue': task_obj.evalue, 'max_hits': task_obj.max_hits},
            outfmt='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'
        )

        if not success_tsv:
            task_obj.status = 'failed'
            task_obj.error_message = blast_error or "BLAST command failed"
            task_obj.save()
            logger.error(f"BLAST TSV analysis for task {task_id} failed: {task_obj.error_message}")
            return JsonResponse({'error': task_obj.error_message}, status=500)

        # 处理原始输出，生成CSV, XML, FASTA等
        try:
            df = pd.read_csv(result_file_path, sep='\t', header=None)
            df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle']
            df.to_csv(csv_result_path, index=False)

            # 运行 -outfmt 5 获取XML报告
            success_xml, xml_error = run_blast_command(
                task_obj.sequence_file,
                requested_database,
                xml_report_path,
                {'program': task_obj.program, 'evalue': task_obj.evalue, 'max_hits': task_obj.max_hits},
                outfmt='5'
            )
            if not success_xml:
                logger.warning(f"Could not generate XML report for task {task_id}. Proceeding without it.")
                task_obj.xml_report_file = None
            else:
                task_obj.xml_report_file = xml_report_path

            # 从TSV结果中提取hit IDs并使用 blastdbcmd 获取对应FASTA序列
            hit_ids = set(df['sseqid'].apply(lambda x: x.split()[0])) # 提取ID部分，去掉描述
            
            # 使用保存的完整数据库路径
            logger.info(f"Fetching FASTA hits using blastdbcmd from {task_obj.database_actual} for task {task_id}")
            fasta_hits = fetch_fasta_hits_via_blastdbcmd(hit_ids, task_obj.database_actual)  # 使用保存的数据库路径
            
            if not fasta_hits:
                 logger.warning(f"No matching sequences found via blastdbcmd for task {task_id} or blastdbcmd failed.")
                 task_obj.fasta_hits_file = None
            else:
                SeqIO.write(fasta_hits, fasta_hits_path, "fasta")
                task_obj.fasta_hits_file = fasta_hits_path # Only store path if hits were found and written

            task_obj.status = 'completed'
            task_obj.result_file = result_file_path
            task_obj.csv_result_file = csv_result_path
            if success_xml:
                task_obj.xml_report_file = xml_report_path
            # fasta_hits_file is set conditionally above
            task_obj.completed_at = timezone.now() # Use Django's timezone-aware now()
            
            task_obj.save() # 保存所有更改到数据库
            
            logger.info(f"BLAST analysis for task {task_id} completed successfully.")
            return JsonResponse({'message': 'Analysis completed successfully'}, status=200)
        except Exception as e:
            logger.error(f"Error processing results for task {task_id}: {str(e)}")
            task_obj.status = 'failed'
            task_obj.error_message = f"Result processing failed: {str(e)}"
            task_obj.save() # 保存失败状态
            return JsonResponse({'error': 'Failed to process results'}, status=500)

    except json.JSONDecodeError:
        logger.error(f"[RUN_ANALYSIS] Invalid JSON in request body for run_analysis on task '{task_id}'. Body: {request.body.decode('utf-8', errors='ignore')}")
        return JsonResponse({'error': 'Invalid JSON in request body for run_analysis'}, status=400)
    except Exception as e:
        logger.error(f"[RUN_ANALYSIS] Error running analysis for task {task_id}: {str(e)}")
        # 尝试更新数据库状态为失败，即使在异常处理中也要尽量更新
        try:
            task_obj.status = 'failed'
            task_obj.error_message = f"Runtime error: {str(e)}"
            task_obj.save()
        except:
            pass # 如果连更新状态都失败了，也没办法了
        return JsonResponse({'error': 'Internal server error during analysis'}, status=500)


def blast_get_task_status(request, task_id):
    if request.method != 'GET':
        return JsonResponse({'error': 'Method not allowed'}, status=405)

    try:
        task_obj = BlastTask.objects.get(id=task_id)
        
        # 使用一个辅助函数安全地格式化日期
        def safe_isoformat(dt):
            if dt:
                try:
                    return dt.isoformat()
                except Exception:
                    # 如果isoformat失败，返回字符串表示
                    return str(dt)
            return None
        
        return JsonResponse({
            'task_id': task_obj.id,
            'status': task_obj.status,
            'name': task_obj.name,
            'database': task_obj.database_requested or 'N/A',
            'created_at': safe_isoformat(task_obj.created_at),
            'completed_at': safe_isoformat(task_obj.completed_at),
            'error_message': task_obj.error_message or '',
        })
        
    except BlastTask.DoesNotExist:
        return JsonResponse({'error': 'Task not found'}, status=404)
    except Exception as e:
        logger.error(f"Error in blast_get_task_status for task {task_id}: {str(e)}")
        return JsonResponse({'error': 'Internal server error'}, status=500)

def blast_get_task_results(request, task_id):
    if request.method != 'GET':
        return JsonResponse({'error': 'Method not allowed'}, status=405)

    # --- 修改：从数据库查询 ---
    try:
        task_obj = BlastTask.objects.get(id=task_id)
    except BlastTask.DoesNotExist:
        return JsonResponse({'error': 'Task not found'}, status=404)
    # --- 修改结束 ---

    if task_obj.status != 'completed':
        return JsonResponse({'error': 'Results not available yet'}, status=400)

    try:
        # 读取 CSV 结果并返回前几行作为预览
        csv_path = task_obj.csv_result_file
        if not csv_path or not os.path.exists(csv_path):
             return JsonResponse({'error': 'Result file not found'}, status=500)

        df = pd.read_csv(csv_path)
        preview_data = df.head(10).to_dict(orient='records') # 返回前10行
        total_hits = len(df)

        # 构建下载链接，只有当文件存在时才提供链接
        download_urls = {
            'csv': f"/hormone_app/api/blast-tasks/{task_id}/download/blast-results/",
        }
        if task_obj.xml_report_file and os.path.exists(task_obj.xml_report_file):
            download_urls['xml'] = f"/hormone_app/api/blast-tasks/{task_id}/download/full-report/"
        if task_obj.fasta_hits_file and os.path.exists(task_obj.fasta_hits_file):
            download_urls['fasta'] = f"/hormone_app/api/blast-tasks/{task_id}/download/fasta-hits/"

        return JsonResponse({
            'task_id': task_id,
            'status': task_obj.status,
            'database': task_obj.database_requested, # 也在这里返回数据库名称
            'total_hits': total_hits,
            'preview': preview_data,
            'download_urls': download_urls, # 只包含存在的文件的下载链接
        })
    except Exception as e:
        logger.error(f"Error retrieving results for task {task_id}: {str(e)}")
        return JsonResponse({'error': 'Failed to retrieve results'}, status=500)


def blast_download_result_csv(request, task_id):
    if request.method != 'GET':
        raise Http404

    # --- 修改：从数据库查询 ---
    try:
        task_obj = BlastTask.objects.get(id=task_id)
    except BlastTask.DoesNotExist:
        raise Http404("Task object not found.")
    # --- 修改结束 ---

    if not task_obj.csv_result_file: # 检查数据库字段
        raise Http404("CSV result file not found.")

    file_path = task_obj.csv_result_file # 从数据库获取路径
    if not os.path.exists(file_path):
        raise Http404("File does not exist on disk.")

    response = FileResponse(
        open(file_path, 'rb'),
        content_type='text/csv'
    )
    response['Content-Disposition'] = f'attachment; filename="{escape_uri_path(os.path.basename(file_path))}"'
    return response

def blast_download_full_report_xml(request, task_id):
    if request.method != 'GET':
        raise Http404

    # --- 修改：从数据库查询 ---
    try:
        task_obj = BlastTask.objects.get(id=task_id)
    except BlastTask.DoesNotExist:
        raise Http404("Task object not found.")
    # --- 修改结束 ---

    # 检查是否存在且文件路径有效
    if not task_obj.xml_report_file or not os.path.exists(task_obj.xml_report_file):
        raise Http404("XML report file not found.")

    file_path = task_obj.xml_report_file
    if not os.path.exists(file_path):
        raise Http404("File does not exist on disk.")

    response = FileResponse(
        open(file_path, 'rb'),
        content_type='application/xml'
    )
    response['Content-Disposition'] = f'attachment; filename="{escape_uri_path(os.path.basename(file_path))}"'
    return response

def blast_download_fasta_hits(request, task_id):
    if request.method != 'GET':
        raise Http404

    # --- 修改：从数据库查询 ---
    try:
        task_obj = BlastTask.objects.get(id=task_id)
    except BlastTask.DoesNotExist:
        raise Http404("Task object not found.")
    # --- 修改结束 ---

    # 检查是否存在且文件路径有效
    if not task_obj.fasta_hits_file or not os.path.exists(task_obj.fasta_hits_file):
        raise Http404("FASTA hits file not found.")

    file_path = task_obj.fasta_hits_file
    if not os.path.exists(file_path):
        raise Http404("File does not exist on disk.")

    response = FileResponse(
        open(file_path, 'rb'),
        content_type='text/plain'
    )
    response['Content-Disposition'] = f'attachment; filename="{escape_uri_path(os.path.basename(file_path))}"'
    return response

def blast_download_all_results(request, task_id):
    if request.method != 'GET':
        raise Http404

    # --- 修改：从数据库查询 ---
    try:
        task_obj = BlastTask.objects.get(id=task_id)
    except BlastTask.DoesNotExist:
        raise Http404("Task object not found.")
    # --- 修改结束 ---

    if task_obj.status != 'completed':
        raise Http404("Task not completed or files not found.")

    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, "a", zipfile.ZIP_DEFLATED, False) as zip_file:
        # 遍历任务信息，只添加存在的文件
        # 注意：这里使用 task_obj 的属性替代原来的 task_info.get(key)
        for attr_name in ['result_file', 'csv_result_file', 'xml_report_file', 'fasta_hits_file']:
            file_attr_value = getattr(task_obj, attr_name, None)
            if file_attr_value and os.path.exists(file_attr_value):
                zip_file.write(file_attr_value, os.path.basename(file_attr_value))
    
    zip_buffer.seek(0)
    response = FileResponse(
        zip_buffer,
        content_type='application/zip'
    )
    response['Content-Disposition'] = f'attachment; filename="blast_task_{task_id}_all_results.zip"'
    return response

def blast_download_evalue_plot(request, task_id):
    if request.method != 'GET':
        raise Http404

    # --- 修改：从数据库查询 ---
    try:
        task_obj = BlastTask.objects.get(id=task_id)
    except BlastTask.DoesNotExist:
        raise Http404("Task object not found.")
    # --- 修改结束 ---

    if task_obj.status != 'completed':
        raise Http404("Task not completed or results not found.")

    try:
        csv_path = task_obj.csv_result_file # 从数据库获取路径
        if not csv_path or not os.path.exists(csv_path):
             raise Http404("CSV result file not found for plotting.")
        
        df = pd.read_csv(csv_path)
        
        plt.figure(figsize=(10, 6))
        plt.hist(df['evalue'], bins=50, log=True, edgecolor='black')
        plt.title('E-value Distribution')
        plt.xlabel('E-value')
        plt.ylabel('Frequency (log scale)')
        plt.grid(True, which="both", ls="-", alpha=0.2)

        buf = BytesIO()
        plt.savefig(buf, format='png')
        plt.close() # 关闭图形以释放内存
        buf.seek(0)

        response = FileResponse(
            buf,
            content_type='image/png'
        )
        response['Content-Disposition'] = f'attachment; filename="evalue_plot_{task_id}.png"'
        return response
    except Exception as e:
        logger.error(f"Error generating evalue plot for task {task_id}: {str(e)}")
        raise Http404("Could not generate plot.")


def blast_download_identity_plot(request, task_id):
    if request.method != 'GET':
        raise Http404

    # --- 修改：从数据库查询 ---
    try:
        task_obj = BlastTask.objects.get(id=task_id)
    except BlastTask.DoesNotExist:
        raise Http404("Task object not found.")
    # --- 修改结束 ---

    if task_obj.status != 'completed':
        raise Http404("Task not completed or results not found.")

    try:
        csv_path = task_obj.csv_result_file # 从数据库获取路径
        if not csv_path or not os.path.exists(csv_path):
             raise Http404("CSV result file not found for plotting.")
        
        df = pd.read_csv(csv_path)
        
        plt.figure(figsize=(10, 6))
        plt.hist(df['pident'], bins=50, edgecolor='black')
        plt.title('Identity Percentage Histogram')
        plt.xlabel('Identity (%)')
        plt.ylabel('Frequency')
        plt.grid(True, which="both", ls="-", alpha=0.2)

        buf = BytesIO()
        plt.savefig(buf, format='png')
        plt.close()
        buf.seek(0)

        response = FileResponse(
            buf,
            content_type='image/png'
        )
        response['Content-Disposition'] = f'attachment; filename="identity_plot_{task_id}.png"'
        return response
    except Exception as e:
        logger.error(f"Error generating identity plot for task {task_id}: {str(e)}")
        raise Http404("Could not generate plot.")


def blast_download_all_visualizations(request, task_id):
    if request.method != 'GET':
        raise Http404

    # --- 修改：从数据库查询 ---
    try:
        task_obj = BlastTask.objects.get(id=task_id)
    except BlastTask.DoesNotExist:
        raise Http404("Task object not found.")
    # --- 修改结束 ---

    if task_obj.status != 'completed':
        raise Http404("Task not completed or results not found.")

    try:
        csv_path = task_obj.csv_result_file # 从数据库获取路径
        if not csv_path or not os.path.exists(csv_path):
             raise Http404("CSV result file not found for plotting.")
        
        df = pd.read_csv(csv_path)

        zip_buffer = BytesIO()

        with zipfile.ZipFile(zip_buffer, "a", zipfile.ZIP_DEFLATED, False) as zip_file:
            
            # Plot 1: E-value
            plt.figure(figsize=(10, 6))
            plt.hist(df['evalue'], bins=50, log=True, edgecolor='black')
            plt.title('E-value Distribution')
            plt.xlabel('E-value')
            plt.ylabel('Frequency (log scale)')
            plt.grid(True, which="both", ls="-", alpha=0.2)
            img_buf = BytesIO()
            plt.savefig(img_buf, format='png')
            plt.close()
            img_buf.seek(0)
            zip_file.writestr(f'evalue_plot_{task_id}.png', img_buf.getvalue())

            # Plot 2: Identity
            plt.figure(figsize=(10, 6))
            plt.hist(df['pident'], bins=50, edgecolor='black')
            plt.title('Identity Percentage Histogram')
            plt.xlabel('Identity (%)')
            plt.ylabel('Frequency')
            plt.grid(True, which="both", ls="-", alpha=0.2)
            img_buf = BytesIO()
            plt.savefig(img_buf, format='png')
            plt.close()
            img_buf.seek(0)
            zip_file.writestr(f'identity_plot_{task_id}.png', img_buf.getvalue())

        zip_buffer.seek(0)
        response = FileResponse(
            zip_buffer,
            content_type='application/zip'
        )
        response['Content-Disposition'] = f'attachment; filename="blast_task_{task_id}_all_visualizations.zip"'
        return response
    except Exception as e:
        logger.error(f"Error generating all visualizations for task {task_id}: {str(e)}")
        raise Http404("Could not generate plots.")
        
def blast_get_alignment_details(request, task_id):
    """获取序列比对详情 - 强制使用正确的查询序列"""
    if request.method != 'GET':
        return JsonResponse({'error': 'Method not allowed'}, status=405)
    
    try:
        task_obj = BlastTask.objects.get(id=task_id)
    except BlastTask.DoesNotExist:
        return JsonResponse({'error': 'Task not found'}, status=404)
    
    if task_obj.status != 'completed':
        return JsonResponse({'error': 'Results not available yet'}, status=400)
    
    try:
        # 1. ✅ 首先获取完整的查询序列（直接从文件读取）
        query_file = task_obj.sequence_file
        if not query_file or not os.path.exists(query_file):
            return JsonResponse({'error': 'Query sequence file not found'}, status=404)
        
        with open(query_file, 'r', encoding='utf-8') as f:
            content = f.read()
            # 提取序列（去掉FASTA头）
            lines = [line.strip() for line in content.split('\n') if not line.startswith('>')]
            full_query_seq = ''.join(lines)
        
        logger.info(f"Query sequence from file: {full_query_seq[:50]}...")
        
        # 2. 读取CSV获取第一个hit的信息
        if not task_obj.csv_result_file or not os.path.exists(task_obj.csv_result_file):
            return JsonResponse({'error': 'CSV result file not found'}, status=404)
        
        df = pd.read_csv(task_obj.csv_result_file)
        if len(df) == 0:
            return JsonResponse({'error': 'No hits found'}, status=404)
        
        first_hit = df.iloc[0]
        
        # 3. 获取比对位置
        qstart = int(first_hit['qstart'])
        qend = int(first_hit['qend'])
        
        # 4. 从查询序列中提取比对片段
        query_fragment = full_query_seq[qstart-1:qend]
        
        # 5. 从XML中获取目标序列的比对片段
        if not task_obj.xml_report_file or not os.path.exists(task_obj.xml_report_file):
            return JsonResponse({'error': 'XML report not found'}, status=404)
        
        with open(task_obj.xml_report_file, 'r', encoding='utf-8') as f:
            xml_content = f.read()
        
        # 提取第一个HSP的hseq
        import re
        hseq_match = re.search(r'<Hsp_hseq>(.*?)</Hsp_hseq>', xml_content, re.DOTALL)
        midline_match = re.search(r'<Hsp_midline>(.*?)</Hsp_midline>', xml_content, re.DOTALL)
        
        if not hseq_match or not midline_match:
            return JsonResponse({'error': 'Could not extract alignment from XML'}, status=404)
        
        hseq = hseq_match.group(1)
        midline = midline_match.group(1)
        
        # 6. 验证长度
        if len(query_fragment) != len(hseq):
            logger.warning(f"Length mismatch: query_fragment={len(query_fragment)}, hseq={len(hseq)}")
            # 截断到较短的长度
            min_len = min(len(query_fragment), len(hseq))
            query_fragment = query_fragment[:min_len]
            hseq = hseq[:min_len]
            midline = midline[:min_len]
        
        # 7. 格式化输出（BLAST经典格式）
        def format_blast_alignment(q, m, s, q_start, width=60):
            lines = []
            
            for i in range(0, len(q), width):
                q_block = q[i:i+width]
                m_block = m[i:i+width]
                s_block = s[i:i+width]
                
                q_pos = q_start + i
                
                lines.append(f"Query  {q_pos:4d}    {q_block}")
                lines.append(f"           {m_block}")
                lines.append(f"Sbjct  {q_pos:4d}    {s_block}")  # 使用相同的position便于阅读
                lines.append("")
            
            return '\n'.join(lines)
        
        alignment_text = format_blast_alignment(query_fragment, midline, hseq, qstart)
        
        # 添加统计信息
        stats = f"""> {first_hit['sseqid']}
Length = {first_hit['length']}
 Score = {first_hit['bitscore']:.1f} bits, Expect = {first_hit['evalue']:.2e}
 Identities = {first_hit['pident']:.1f}% ({int(first_hit['length'] * first_hit['pident']/100)}/{first_hit['length']})
 
"""
        
        # 添加调试信息
        debug_info = f"""
Debug Info:
- Query file: {query_file}
- Query sequence length: {len(full_query_seq)}
- Query fragment: {qstart}-{qend} ({len(query_fragment)} aa)
- First 20 aa of query: {full_query_seq[:20]}
- First 20 aa of fragment: {query_fragment[:20]}
"""
        
        return JsonResponse({
            'alignment': stats + alignment_text,
            'debug': debug_info,
            'task_id': task_id
        })
        
    except Exception as e:
        logger.error(f"Error generating alignment: {str(e)}")
        return JsonResponse({'error': str(e)}, status=500)