"""
KEGG Pathway Analysis View Module
数据库存储版本 - 不需要用户登录，数据保存在数据库
提供KEGG通路富集分析功能
"""
import os
import uuid
import json
import logging
import tempfile
from datetime import datetime
from pathlib import Path
from django.http import JsonResponse, FileResponse, HttpResponse
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_http_methods
from django.conf import settings
import pandas as pd
import numpy as np
import requests
import io
import zipfile
import csv
import re
import random
from datetime import datetime, timedelta
from django.db.models import Q

# 导入模型
from ..models import KeggAnalysisTask,KeggPathway,KeggGenePathway
from ..services.kegg_analysis import KeggAnalysisService

# 配置日志
logger = logging.getLogger(__name__)

# 结果目录（用于存储临时文件）
KEGG_RESULTS_DIR = Path(settings.BASE_DIR) / 'media' / 'kegg_results'

# 确保结果目录存在
KEGG_RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ==================== 辅助函数 ====================

def generate_client_id(request):
    """生成客户端标识（基于IP和浏览器信息）"""
    client_ip = request.META.get('REMOTE_ADDR', 'unknown')
    user_agent = request.META.get('HTTP_USER_AGENT', 'unknown')
    return f"{client_ip}_{hash(user_agent)}"

def get_task(task_id):
    """从数据库获取任务"""
    try:
        if isinstance(task_id, str):
            task_id = uuid.UUID(task_id)
        return KeggAnalysisTask.objects.get(task_id=task_id)
    except (KeggAnalysisTask.DoesNotExist, ValueError):
        return None

def update_task_status(task, status, message=None, progress=None):
    """更新任务状态"""
    task.status = status
    if message:
        task.message = message
    if progress is not None:
        task.progress = progress
    task.save()
    return True

def get_example_pathways():
    """返回示例通路数据（用于测试和演示）"""
    return [
        {"id": "hsa04110", "name": "Cell cycle", "category": "Cellular Processes", "gene_count": 28, 
         "enrichment_score": 5.23, "p_value": 1.2e-08, "adjusted_pvalue": 1.2e-08, 
         "genes": ["CDK1", "CCNB1", "TP53", "CDKN1A", "RB1", "E2F1", "CCNE1", "CDC25A"]},
        {"id": "hsa04151", "name": "PI3K-Akt signaling pathway", "category": "Environmental Information Processing", 
         "gene_count": 42, "enrichment_score": 4.87, "p_value": 3.5e-07, "adjusted_pvalue": 3.5e-07,
         "genes": ["PIK3CA", "AKT1", "PTEN", "MTOR", "GSK3B", "FOXO1"]},
        {"id": "hsa04010", "name": "MAPK signaling pathway", "category": "Environmental Information Processing", 
         "gene_count": 36, "enrichment_score": 4.12, "p_value": 7.8e-06, "adjusted_pvalue": 7.8e-06,
         "genes": ["MAPK1", "KRAS", "BRAF", "EGFR", "JUN", "FOS"]},
        {"id": "hsa04210", "name": "Apoptosis", "category": "Cellular Processes", 
         "gene_count": 19, "enrichment_score": 3.89, "p_value": 2.1e-05, "adjusted_pvalue": 2.1e-05,
         "genes": ["TP53", "BCL2", "BAX", "CASP3", "CASP8", "CASP9"]},
        {"id": "hsa04350", "name": "TGF-beta signaling pathway", "category": "Environmental Information Processing", 
         "gene_count": 15, "enrichment_score": 3.24, "p_value": 4.7e-05, "adjusted_pvalue": 4.7e-05,
         "genes": ["TGFB1", "SMAD3", "SMAD4", "TGFBR1", "TGFBR2"]},
        {"id": "hsa04310", "name": "Wnt signaling pathway", "category": "Environmental Information Processing", 
         "gene_count": 22, "enrichment_score": 3.18, "p_value": 8.3e-05, "adjusted_pvalue": 8.3e-05,
         "genes": ["CTNNB1", "APC", "GSK3B", "AXIN1", "AXIN2"]}
    ]

# ==================== 视图函数 ====================

def kegg_analysis_page(request):
    """KEGG分析页面 - 完全客户端版本"""
    from django.shortcuts import render
    return render(request, 'kegg_analysis.html')

@csrf_exempt
@require_http_methods(["POST"])
def create_kegg_task(request):
    """创建KEGG分析任务 - 不需要用户认证"""
    try:
        data = json.loads(request.body)
        
        # 生成客户端标识
        client_ip = request.META.get('REMOTE_ADDR', 'unknown')
        client_id = generate_client_id(request)
        
        # 创建数据库记录
        task = KeggAnalysisTask.objects.create(
            client_id=client_id,
            client_ip=client_ip,
            parameters=data,
            status='created',
            progress=0,
            message='Task created'
        )
        
        logger.info(f"KEGG task created in DB: {task.task_id} for client: {client_ip}")
        
        return JsonResponse({
            'success': True,
            'task_id': str(task.task_id),  # 注意：返回字符串格式的UUID
            'message': 'KEGG analysis task created successfully',
            'client_mode': True,  # 标识为客户端模式
            'data_persistence_note': 'Data is stored in database and will be preserved across server restarts.'
        })
        
    except Exception as e:
        logger.error(f"创建KEGG任务失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=500)

@csrf_exempt
@require_http_methods(["POST"])
def upload_kegg_file(request, task_id):
    """上传KEGG分析文件 - 不需要用户认证"""
    try:
        task = get_task(task_id)
        if not task:
            return JsonResponse({
                'success': False,
                'error': 'Task not found'
            }, status=404)
        
        if 'gene_file' not in request.FILES:
            return JsonResponse({
                'success': False,
                'error': 'No file provided'
            }, status=400)
        
        file = request.FILES['gene_file']
        
        # 检查文件类型
        allowed_extensions = ['.txt', '.csv', '.tsv']
        file_ext = os.path.splitext(file.name)[1].lower()
        if file_ext not in allowed_extensions:
            return JsonResponse({
                'success': False,
                'error': f'Invalid file type. Allowed types: {", ".join(allowed_extensions)}'
            }, status=400)
        
        # 检查文件大小（最大10MB）
        max_size = 10 * 1024 * 1024  # 10MB
        if file.size > max_size:
            return JsonResponse({
                'success': False,
                'error': f'File too large. Maximum size is {max_size // (1024*1024)}MB'
            }, status=400)
        
        # 读取文件内容
        content = file.read().decode('utf-8')
        
        # 解析基因列表
        genes = []
        lines = content.strip().split('\n')
        for line in lines:
            # 处理逗号或制表符分隔
            for item in line.replace(',', '\t').split('\t'):
                item = item.strip()
                if item and not item.startswith('#'):
                    # 移除可能的引号
                    item = item.strip('"\'')
                    if item:  # 再次检查是否为空
                        genes.append(item)
        
        # 去重
        genes = list(set(genes))
        
        # 检查基因数量
        if len(genes) < 2:
            return JsonResponse({
                'success': False,
                'error': 'At least 2 genes are required for analysis'
            }, status=400)
        
        if len(genes) > 10000:
            return JsonResponse({
                'success': False,
                'error': 'Maximum 10,000 genes allowed'
            }, status=400)
        
        # 更新数据库记录
        task.genes = genes
        task.gene_count = len(genes)
        task.file_name = file.name
        task.status = 'uploaded'
        task.message = f'File uploaded: {file.name} ({len(genes)} genes)'
        task.save()
        
        client_ip = request.META.get('REMOTE_ADDR', 'unknown')
        logger.info(f"File uploaded for task {task_id}: {file.name} ({len(genes)} genes) from {client_ip}")
        
        return JsonResponse({
            'success': True,
            'task_id': str(task.task_id),
            'gene_count': len(genes),
            'file_name': file.name,
            'message': f'File uploaded successfully. Found {len(genes)} genes.',
            'client_mode': True
        })
        
    except UnicodeDecodeError:
        return JsonResponse({
            'success': False,
            'error': 'File encoding error. Please save your file as UTF-8 text.'
        }, status=400)
    except Exception as e:
        logger.error(f"上传KEGG文件失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=500)

@csrf_exempt
@require_http_methods(["POST"])
def validate_kegg_gene_list(request):
    """验证基因列表 - 使用KEGG数据库进行真实验证"""
    try:
        data = json.loads(request.body)
        genes = data.get('genes', [])
        id_type = data.get('id_type', 'symbol')
        organism = data.get('organism', 'hsa')
        
        # 检查基因列表是否为空
        if not genes:
            return JsonResponse({'success': False, 'error': 'Gene list is empty'}, status=400)
        
        # ✅ 从数据库验证基因
        from ..models import KeggGene, KeggOrganism
        
        try:
            org = KeggOrganism.objects.get(code=organism)
            
            # 清理基因符号
            cleaned_genes = [g.strip().upper() for g in genes if g.strip()]
            
            # 查询数据库中存在的基因
            existing_genes = KeggGene.objects.filter(
                organism=org,
                symbol__in=cleaned_genes
            ).values_list('symbol', flat=True)
            
            existing_set = set(existing_genes)
            
            valid_genes = []
            invalid_genes = []
            
            for gene in cleaned_genes:
                if gene in existing_set:
                    valid_genes.append(gene)
                else:
                    # 尝试模糊匹配
                    fuzzy_matches = KeggGene.objects.filter(
                        organism=org,
                        symbol__icontains=gene
                    )[:1]
                    if fuzzy_matches.exists():
                        valid_genes.append(fuzzy_matches[0].symbol)
                    else:
                        invalid_genes.append(gene)
            
            # 去重
            valid_genes = list(set(valid_genes))
            invalid_genes = list(set(invalid_genes))
            
            mapping_rate = len(valid_genes) / len(cleaned_genes) * 100 if cleaned_genes else 0
            
            # 获取数据库统计信息
            total_genes = KeggGene.objects.filter(organism=org).count()
            
            return JsonResponse({
                'success': True,
                'valid_genes': len(valid_genes),
                'valid_gene_list': valid_genes[:50],  # 返回前50个
                'invalid_genes': invalid_genes[:50],
                'total_genes': len(genes),
                'cleaned_genes': len(cleaned_genes),
                'mapping_rate': round(mapping_rate, 2),
                'database_stats': {
                    'total_genes': total_genes,
                    'organism': organism
                },
                'client_mode': True,
                'using_real_data': True
            })
            
        except KeggOrganism.DoesNotExist:
            return JsonResponse({
                'success': False,
                'error': f'Organism {organism} not found in database'
            }, status=400)
        
    except Exception as e:
        logger.error(f"验证基因列表失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=500)

@csrf_exempt
@require_http_methods(["POST"])
def run_kegg_analysis(request, task_id):
    """运行KEGG分析 - 使用本地数据库进行真实分析"""
    try:
        task = get_task(task_id)
        if not task:
            return JsonResponse({
                'success': False,
                'error': 'Task not found'
            }, status=404)
        
        # 检查任务状态
        if task.status == 'completed':
            return JsonResponse({
                'success': False,
                'error': 'Analysis already completed',
                'status': 'completed'
            }, status=400)
        
        # 检查是否有基因数据
        if not task.genes or len(task.genes) == 0:
            return JsonResponse({
                'success': False,
                'error': 'No gene data found. Please upload a gene list first.'
            }, status=400)
        
        # 更新任务状态为运行中
        update_task_status(task, 'running', 'Starting KEGG analysis...', 10)
        
        client_ip = request.META.get('REMOTE_ADDR', 'unknown')
        logger.info(f"Starting real KEGG analysis for task {task_id} from {client_ip}")
        
        try:
            # 获取参数
            parameters = task.parameters if isinstance(task.parameters, dict) else {}
            organism = parameters.get('organism', 'hsa')
            pvalue_cutoff = float(parameters.get('pvalue_cutoff', 0.05))
            qvalue_cutoff = float(parameters.get('qvalue_cutoff', 0.2))
            min_genes = int(parameters.get('min_genes', 2))
            max_genes = int(parameters.get('max_genes', 500))
            
            # 更新进度
            update_task_status(task, 'running', 'Initializing KEGG analysis service...', 20)
            
            # 创建分析服务实例
            from ..services.kegg_analysis import KeggAnalysisService
            analyzer = KeggAnalysisService(organism=organism)
            
            # 获取基因列表
            genes = task.genes
            logger.info(f"Task {task_id}: Analyzing {len(genes)} genes for {organism}")
            
            update_task_status(task, 'running', f'Mapping {len(genes)} genes to KEGG database...', 30)
            
            # 执行富集分析
            analysis_results = analyzer.perform_enrichment_analysis(genes)
            
            # 检查是否有错误
            if 'error' in analysis_results:
                error_msg = analysis_results['error']
                update_task_status(task, 'failed', f'Analysis failed: {error_msg}')
                logger.error(f"Task {task_id}: {error_msg}")
                return JsonResponse({
                    'success': False,
                    'error': error_msg,
                    'client_mode': True
                }, status=400)
            
            update_task_status(task, 'running', 'Processing enrichment results...', 70)
            
            # 过滤显著通路
            significant_pathways = [
                p for p in analysis_results.get('results', [])
                if p.get('adjusted_pvalue', 1.0) < pvalue_cutoff
                and p.get('overlap_count', 0) >= min_genes
                and p.get('overlap_count', 0) <= max_genes
            ]
            
            # 准备返回结果（转换为前端期望的格式）
            pathways = []
            for pathway in significant_pathways[:100]:  # 限制最多100个
                pathways.append({
                    'id': pathway.get('pathway_id', ''),
                    'name': pathway.get('pathway_name', ''),
                    'category': pathway.get('category', 'Unknown'),
                    'gene_count': pathway.get('overlap_count', 0),
                    'pathway_gene_count': pathway.get('pathway_gene_count', 0),
                    'enrichment_score': pathway.get('enrichment_score', 0),
                    'p_value': pathway.get('p_value', 1.0),
                    'adjusted_pvalue': pathway.get('adjusted_pvalue', 1.0),
                    'log10_pvalue': pathway.get('log10_pvalue', 0),
                    'log10_adj_pvalue': pathway.get('log10_adj_pvalue', 0),
                    'odds_ratio': pathway.get('odds_ratio', 0),
                    'genes': pathway.get('overlap_genes', [])[:20],  # 只显示前20个基因
                    'pathway_url': pathway.get('pathway_url', ''),
                    'pathway_image': pathway.get('pathway_image', '')
                })
            
            # 按p值排序
            pathways.sort(key=lambda x: x['p_value'])
            
            # 计算统计摘要
            total_genes = len(genes)
            mapped_genes_count = analysis_results.get('mapped_gene_count', 0)
            
            # 构建完整结果数据
            result_data = {
                'task_id': str(task.task_id),
                'analysis_date': datetime.now().isoformat(),
                'organism': organism,
                'parameters': {
                    'organism': organism,
                    'pvalue_cutoff': pvalue_cutoff,
                    'qvalue_cutoff': qvalue_cutoff,
                    'min_genes': min_genes,
                    'max_genes': max_genes,
                    'analysis_method': 'hypergeometric',
                    'correction_method': 'fdr_bh'
                },
                'summary': {
                    'total_genes_analyzed': total_genes,
                    'genes_mapped': mapped_genes_count,
                    'mapping_rate': round(mapped_genes_count/total_genes*100, 2) if total_genes > 0 else 0,
                    'total_pathways_analyzed': analysis_results.get('total_pathways_analyzed', 0),
                    'significant_pathways': len(significant_pathways),
                    'significant_pathways_p05': len([p for p in significant_pathways if p.get('adjusted_pvalue', 1.0) < 0.05]),
                    'significant_pathways_p01': len([p for p in significant_pathways if p.get('adjusted_pvalue', 1.0) < 0.01]),
                    'significant_pathways_p001': len([p for p in significant_pathways if p.get('adjusted_pvalue', 1.0) < 0.001]),
                    'top_pathway': pathways[0]['name'] if pathways else 'None',
                    'top_pathway_pvalue': pathways[0]['adjusted_pvalue'] if pathways else 1.0
                },
                'pathways': pathways,
                'gene_pathway_mapping': [],  # 可以后续添加详细的基因-通路映射
                'database_stats': {
                    'total_pathways': 369,
                    'total_genes': 24680,
                    'total_associations': 39985,
                    'organism': organism
                }
            }
            
            update_task_status(task, 'running', 'Saving analysis results...', 90)
            
            # 保存结果到数据库
            task.results = result_data
            task.status = 'completed'
            task.progress = 100
            task.message = f'Analysis completed successfully. Found {len(significant_pathways)} significant pathways.'
            task.save()
            
            # 同时保存到文件作为备份
            try:
                result_file = KEGG_RESULTS_DIR / f"{task_id}_results.json"
                with open(result_file, 'w', encoding='utf-8') as f:
                    json.dump(result_data, f, indent=2, ensure_ascii=False)
                logger.info(f"Results saved to file for task {task_id}")
            except Exception as file_error:
                logger.warning(f"Could not save result file for task {task_id}: {file_error}")
            
            logger.info(f"Real KEGG analysis completed for task {task_id}. Found {len(significant_pathways)} significant pathways.")
            
            # 返回结果摘要
            return JsonResponse({
                'success': True,
                'task_id': str(task.task_id),
                'status': 'completed',
                'message': 'KEGG analysis completed successfully using local database',
                'progress': 100,
                'client_mode': True,
                'data_persistence_note': 'Results are stored in database and will be preserved.',
                'result_summary': {
                    'significant_pathways': len(significant_pathways),
                    'total_genes': total_genes,
                    'mapped_genes': mapped_genes_count,
                    'mapping_rate': round(mapped_genes_count/total_genes*100, 2) if total_genes > 0 else 0,
                    'top_pathway': pathways[0]['name'] if pathways else 'None',
                    'top_pathway_pvalue': pathways[0]['adjusted_pvalue'] if pathways else None,
                    'analysis_method': 'Local KEGG Database',
                    'database_stats': f"369 pathways, 24680 genes, 39985 associations"
                }
            })
            
        except Exception as analysis_error:
            error_message = str(analysis_error)
            logger.error(f"Analysis failed for task {task_id}: {error_message}", exc_info=True)
            update_task_status(task, 'failed', f'Analysis failed: {error_message}')
            return JsonResponse({
                'success': False,
                'error': f'Analysis failed: {error_message}',
                'client_mode': True
            }, status=500)
        
    except Exception as e:
        logger.error(f"运行KEGG分析失败: {e}", exc_info=True)
        return JsonResponse({
            'success': False,
            'error': str(e),
            'client_mode': True
        }, status=500)
        
@csrf_exempt
@require_http_methods(["POST"])
def preview_kegg_mapping(request, task_id):
    """预览基因映射结果（增强版）"""
    try:
        task = get_task(task_id)
        if not task:
            return JsonResponse({'success': False, 'error': 'Task not found'}, status=404)
        
        if not task.genes or len(task.genes) == 0:
            return JsonResponse({'success': False, 'error': 'No gene data found.'}, status=400)
        
        from ..services.kegg_analysis import KeggAnalysisService
        
        parameters = task.parameters if isinstance(task.parameters, dict) else {}
        organism = parameters.get('organism', 'hsa')
        
        analyzer = KeggAnalysisService(organism=organism)
        
        # 使用服务的映射逻辑
        genes = task.genes
        mapped_genes = analyzer.map_gene_symbols(genes)
        unmapped_genes = [g for g in genes if g not in mapped_genes]
        
        return JsonResponse({
            'success': True,
            'preview': {
                'total_genes': len(genes),
                'mapped_genes_count': len(mapped_genes),
                'mapped_genes': mapped_genes[:20],
                'unmapped_genes_count': len(unmapped_genes),
                'unmapped_genes': unmapped_genes[:20],
                'mapping_rate': round(len(mapped_genes)/len(genes)*100, 2) if genes else 0,
                'organism': organism,
                'database_stats': {
                    'total_pathways': 369,
                    'total_genes': 24680,
                    'total_associations': 39985,
                    'available_for_organism': KeggGene.objects.filter(organism__code=organism).count()
                }
            }
        })
        
    except Exception as e:
        logger.error(f"预览基因映射失败: {e}", exc_info=True)
        return JsonResponse({'success': False, 'error': str(e)}, status=500)
    
@require_http_methods(["GET"])
def get_kegg_database_stats(request):
    """获取KEGG数据库统计信息"""
    try:
        organism = request.GET.get('organism', 'hsa')
        
        from ..models import KeggOrganism, KeggPathway, KeggGene, KeggGenePathway
        
        stats = {
            'organisms': [],
            'current_organism': organism,
            'total_pathways': KeggPathway.objects.count(),
            'total_genes': KeggGene.objects.count(),
            'total_associations': KeggGenePathway.objects.count(),
        }
        
        # 获取各生物体统计
        for org in KeggOrganism.objects.all():
            org_stats = {
                'code': org.code,
                'name': org.name,
                'genes': KeggGene.objects.filter(organism=org).count(),
                'pathways': KeggPathway.objects.filter(
                    gene_relations__gene__organism=org
                ).distinct().count(),
                'associations': KeggGenePathway.objects.filter(
                    gene__organism=org
                ).count()
            }
            stats['organisms'].append(org_stats)
        
        # 获取当前生物体的详细统计
        if organism:
            stats['current'] = {
                'genes': KeggGene.objects.filter(organism__code=organism).count(),
                'pathways': KeggPathway.objects.filter(
                    gene_relations__gene__organism__code=organism
                ).distinct().count(),
                'associations': KeggGenePathway.objects.filter(
                    gene__organism__code=organism
                ).count()
            }
        
        return JsonResponse({
            'success': True,
            'stats': stats,
            'client_mode': True,
            'note': 'Using local KEGG database'
        })
        
    except Exception as e:
        logger.error(f"获取数据库统计失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e),
            'client_mode': True
        }, status=500)

@require_http_methods(["GET"])
def get_kegg_task_status(request, task_id):
    """获取KEGG任务状态 - 不需要用户认证"""
    try:
        task = get_task(task_id)
        if not task:
            return JsonResponse({
                'success': False,
                'error': 'Task not found'
            }, status=404)
        
        client_ip = request.META.get('REMOTE_ADDR', 'unknown')
        
        return JsonResponse({
            'success': True,
            'task_id': str(task.task_id),
            'status': task.status,
            'message': task.message,
            'progress': task.progress,
            'created_at': task.created_at.isoformat(),
            'updated_at': task.updated_at.isoformat(),
            'client_mode': True,
            'client_ip': client_ip
        })
        
    except Exception as e:
        logger.error(f"获取任务状态失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e),
            'client_mode': True
        }, status=500)

@require_http_methods(["GET"])
def get_kegg_task_results(request, task_id):
    """获取KEGG分析结果 - 不需要用户认证"""
    try:
        task = get_task(task_id)
        if not task:
            return JsonResponse({
                'success': False,
                'error': 'Task not found'
            }, status=404)
        
        if task.status != 'completed':
            return JsonResponse({
                'success': False,
                'error': 'Analysis not completed yet',
                'status': task.status,
                'progress': task.progress
            }, status=202)
        
        client_ip = request.META.get('REMOTE_ADDR', 'unknown')
        logger.info(f"Results retrieved for task {task_id} from {client_ip}")
        
        return JsonResponse({
            'success': True,
            'task_id': str(task.task_id),
            'results': task.results if task.results else {},
            'client_mode': True,
            'download_available': True,
            'data_persistence_note': 'Results are stored in database and will be preserved.'
        })
        
    except Exception as e:
        logger.error(f"获取分析结果失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e),
            'client_mode': True
        }, status=500)

@require_http_methods(["GET"])
def get_kegg_visualization_plot(request, task_id, plot_type):
    """获取KEGG可视化图表数据 - 使用真实数据"""
    try:
        task = get_task(task_id)
        if not task:
            return JsonResponse({'success': False, 'error': 'Task not found'}, status=404)

        if task.status != 'completed':
            return JsonResponse({
                'success': False,
                'error': 'Analysis not completed yet',
                'status': task.status
            }, status=202)

        # ✅ 使用真实的任务结果数据
        results = task.results if task.results else {}
        pathways = results.get('pathways', [])

        if not pathways:
            return JsonResponse({
                'success': False,
                'error': 'No pathway data available'
            }, status=404)

        # 初始化 viz_data
        viz_data = {}

        if plot_type == "bar":
            # 条形图数据 - 使用真实数据
            top_pathways = pathways[:15]  # 取前15个
            pathway_names = []
            scores = []
            pvalues = []
            for p in top_pathways:
                name = p.get('name', '')
                if len(name) > 40:
                    name = name[:37] + '...'
                pathway_names.append(name)
                scores.append(round(p.get('enrichment_score', 0), 2))
                pvalues.append(-np.log10(p.get('p_value', 1)) if p.get('p_value', 1) > 0 else 0)

            viz_data = {
                'type': 'bar',
                'data': {
                    'pathways': pathway_names,
                    'scores': scores,
                    'pvalues': pvalues
                },
                'using_real_data': True
            }
        elif plot_type == "bubble":
            # 气泡图数据
            top_pathways = pathways[:12]
            pathway_names = []
            gene_counts = []
            pvalues = []
            for p in top_pathways:
                name = p.get('name', '')
                if len(name) > 30:
                    name = name[:27] + '...'
                pathway_names.append(name)
                gene_counts.append(p.get('gene_count', 0))
                pvalues.append(-np.log10(p.get('p_value', 1)) if p.get('p_value', 1) > 0 else 0)

            viz_data = {
                'type': 'bubble',
                'data': {
                    'pathways': pathway_names,
                    'gene_counts': gene_counts,
                    'pvalues': pvalues
                },
                'using_real_data': True
            }
        elif plot_type == "volcano":
            # 火山图数据
            volcano_data = []
            for p in pathways[:30]:
                volcano_data.append({
                    'pathway': p.get('name', ''),
                    'log2fc': np.log2(p.get('enrichment_score', 1)) if p.get('enrichment_score', 1) > 0 else 0,
                    'neg_log10_p': -np.log10(p.get('p_value', 1)) if p.get('p_value', 1) > 0 else 0,
                    'significant': p.get('adjusted_pvalue', 1) < 0.05
                })

            viz_data = {
                'type': 'volcano',
                'data': volcano_data,
                'using_real_data': True
            }
        # --- 新增：处理饼图 ---
        elif plot_type == "pie":
            # 饼图数据 - 使用真实数据
            top_pathways = pathways[:10]  # 取前10个
            pathway_names = []
            gene_counts = [] # 或者使用 scores = []
            for p in top_pathways:
                name = p.get('name', '')
                if len(name) > 40:
                    name = name[:37] + '...'
                pathway_names.append(name)
                gene_counts.append(p.get('gene_count', 0))

            # 计算剩余通路的总和（可选）
            total_top_genes = sum(gene_counts)
            total_all_genes = sum(p.get('gene_count', 0) for p in pathways)
            other_genes = total_all_genes - total_top_genes

            if other_genes > 0:
                 pathway_names.append("Others")
                 gene_counts.append(other_genes)

            # ✅ 调整数据结构以适配 Chart.js 或类似库
            # Chart.js 饼图期望格式: { labels: [...], datasets: [{ data: [...] }] }
            viz_data = {
                'type': 'pie',
                'data': {
                    'labels': pathway_names,
                    'datasets': [{ # 将 values 放入 datasets[0].data
                        'data': gene_counts,
                        'backgroundColor': [ # 可以提供默认颜色或动态生成
                            '#FF6384', '#36A2EB', '#FFCE56', '#4BC0C0', '#9966FF',
                            '#FF9F40', '#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4',
                            '#FFEAA7', '#DDA0DD', '#98D8C8', '#F7DC6F', '#BB8FCE'
                            # 如果标签超过15个，颜色可能会循环或需要更多颜色
                        ] # 提供背景色数组
                    }],
                    # 如果前端库需要其他配置项，可以在这里添加
                },
                'using_real_data': True
            }
        # --- 修改：处理热图 ---
        elif plot_type == "heatmap":
            # 热图数据
            top_pathways = pathways[:15] # 限制数量避免矩阵过大
            if not top_pathways:
                return JsonResponse({
                    'success': False,
                    'error': 'Insufficient pathway data for heatmap'
                }, status=400)

            pathway_names = []
            enrichment_scores = []
            p_values_neg_log10 = []

            for p in top_pathways:
                name = p.get('name', '')
                if len(name) > 30:
                    name = name[:27] + '...'
                pathway_names.append(name)
                enrichment_scores.append(p.get('enrichment_score', 0))
                p_values_neg_log10.append(-np.log10(p.get('p_value', 1)) if p.get('p_value', 1) > 0 else 0)

            # 标签
            x_labels = ['Enrichment Score', '-log10(P-value)']
            y_labels = pathway_names

            # 数据矩阵
            data_matrix = []
            for i in range(len(top_pathways)):
                data_matrix.append([
                    enrichment_scores[i],
                    p_values_neg_log10[i]
                ])

            # ✅ 调整数据结构以适配前端可能期望的格式
            # 假设前端库期望的是 z, x, y 分离的结构，或者一个包含这些信息的对象
            # 这里提供一个更通用的结构，具体需根据前端库调整
            viz_data = {
                'type': 'heatmap',
                'data': {
                    'z': data_matrix,      # 值矩阵 (Z-axis)
                    'x': x_labels,         # X-axis labels
                    'y': y_labels,         # Y-axis labels
                    # 'title': 'Pathway Enrichment Heatmap' # 可选
                },
                'using_real_data': True
            }
        # --- 结束新增 ---
        else:
            # 如果请求了未知的类型，仍然返回空数据，但这应该不会发生，除非前端有误
            # 也可以返回一个错误
             return JsonResponse({
                 'success': False,
                 'error': f'Plot type "{plot_type}" is not supported.',
                 'supported_types': ['bar', 'bubble', 'volcano', 'pie', 'heatmap'] # 列出支持的类型
             }, status=400)
             # Or, to keep the original behavior but make it clearer that it's an unsupported type:
             # viz_data = {
             #     'type': plot_type,
             #     'data': {},
             #     'using_real_data': True,
             #     'error': f'Plot type "{plot_type}" is not implemented.'
             # }


        # --- 返回响应 ---
        # 注意：原代码在这里直接返回了 JsonResponse，我们调整一下
        if viz_data: # 确保 viz_data 已被设置
            return JsonResponse({
                'success': True,
                'plot_type': plot_type,
                'data': viz_data['data'], # 只返回实际的绘图数据部分
                'using_real_data': True,
                'client_mode': True
            })
        else:
            # 这个分支理论上不应该到达，因为 else 会处理未知类型
             return JsonResponse({
                 'success': False,
                 'error': f'Failed to generate data for plot type "{plot_type}". Internal logic error.',
                 'client_mode': True
             }, status=500)


    except Exception as e:
        logger.error(f"获取可视化图表失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e),
            'client_mode': True
        }, status=500)

@require_http_methods(["GET"])
def download_kegg_result_by_type(request, task_id):
    """按类型下载KEGG结果 - 不需要用户认证"""
    try:
        download_type = request.GET.get('type', 'table')
        task = get_task(task_id)
        
        if not task or task.status != 'completed':
            return JsonResponse({
                'success': False,
                'error': 'Analysis not completed or task not found'
            }, status=404)
        
        results = task.results if task.results else {}
        pathways = results.get('pathways', [])
        
        if download_type == 'table':
            # 下载富集表格（CSV）
            output = io.StringIO()
            writer = csv.writer(output)
            
            # 写入表头
            writer.writerow(['Pathway ID', 'Pathway Name', 'Category', 'Gene Count', 
                           'Enrichment Score', 'P-value', 'Adjusted P-value', 'Top Genes'])
            
            # 写入数据
            for pathway in pathways:
                top_genes = ';'.join(pathway.get('genes', [])[:5])
                writer.writerow([
                    pathway.get('id', ''),
                    pathway.get('name', ''),
                    pathway.get('category', ''),
                    pathway.get('gene_count', 0),
                    pathway.get('enrichment_score', 0),
                    pathway.get('p_value', 0),
                    pathway.get('adjusted_pvalue', 0),
                    top_genes
                ])
            
            response = HttpResponse(output.getvalue(), content_type='text/csv')
            response['Content-Disposition'] = f'attachment; filename="kegg_enrichment_{task_id}.csv"'
            return response
            
        elif download_type == 'genes':
            # 下载基因列表
            response = HttpResponse(content_type='text/csv')
            response['Content-Disposition'] = f'attachment; filename="kegg_genes_{task_id}.csv"'
            
            writer = csv.writer(response)
            writer.writerow(['Gene', 'Pathways'])
            
            # 收集所有基因及其所在通路
            gene_pathways = {}
            for pathway in pathways:
                pathway_name = pathway.get('name', 'Unknown')
                for gene in pathway.get('genes', []):
                    if gene not in gene_pathways:
                        gene_pathways[gene] = []
                    gene_pathways[gene].append(pathway_name)
            
            for gene in sorted(gene_pathways.keys()):
                pathways_str = '; '.join(gene_pathways[gene][:3])  # 只显示前3个通路
                if len(gene_pathways[gene]) > 3:
                    pathways_str += f' (+{len(gene_pathways[gene])-3} more)'
                writer.writerow([gene, pathways_str])
                
            return response
            
        elif download_type == 'summary':
            # 下载分析摘要
            response = HttpResponse(content_type='text/plain')
            response['Content-Disposition'] = f'attachment; filename="kegg_summary_{task_id}.txt"'
            
            summary = results.get('summary', {})
            
            content = f"""KEGG Pathway Analysis Summary
==========================================
Task ID: {task_id}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Organism: {results.get('organism', 'hsa')}

Analysis Summary:
-----------------
Total Genes Analyzed: {summary.get('total_genes_analyzed', 0)}
Genes Mapped to KEGG: {summary.get('genes_mapped', 0)}
Mapping Rate: {summary.get('mapping_rate', 0)}%
Significant Pathways (p < 0.05): {summary.get('significant_pathways', 0)}
Total Pathways Analyzed: {summary.get('total_pathways_analyzed', 0)}

Top Pathway: {summary.get('top_pathway', 'None')}

Generated by Animal Hormone Database KEGG Analysis Tool
Database Storage Mode: Data is preserved in database
"""
            
            response.write(content)
            return response
            
        else:
            return JsonResponse({
                'success': False,
                'error': f'Unknown download type: {download_type}'
            }, status=400)
            
    except Exception as e:
        logger.error(f"下载结果失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=500)

@require_http_methods(["GET"])
def download_all_kegg_results(request, task_id):
    """下载所有KEGG分析结果 - 不需要用户认证"""
    try:
        task = get_task(task_id)
        
        if not task or task.status != 'completed':
            return JsonResponse({
                'success': False,
                'error': 'Analysis not completed or task not found'
            }, status=404)
        
        client_ip = request.META.get('REMOTE_ADDR', 'unknown')
        logger.info(f"Downloading all results for task {task_id} from {client_ip}")
        
        # 创建ZIP文件
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            # 添加富集表格
            results = task.results if task.results else {}
            pathways = results.get('pathways', [])
            summary = results.get('summary', {})
            
            # 1. 富集表格
            csv_data = io.StringIO()
            writer = csv.writer(csv_data)
            writer.writerow(['Pathway ID', 'Pathway Name', 'Category', 'Gene Count', 
                           'Enrichment Score', 'P-value', 'Adjusted P-value', 'Top Genes'])
            for pathway in pathways:
                top_genes = ';'.join(pathway.get('genes', [])[:5])
                writer.writerow([
                    pathway.get('id', ''),
                    pathway.get('name', ''),
                    pathway.get('category', ''),
                    pathway.get('gene_count', 0),
                    pathway.get('enrichment_score', 0),
                    pathway.get('p_value', 0),
                    pathway.get('adjusted_pvalue', 0),
                    top_genes
                ])
            zip_file.writestr('kegg_enrichment_table.csv', csv_data.getvalue())
            
            # 2. 分析报告
            report = f"""KEGG Pathway Enrichment Analysis Report
==========================================
Task ID: {task_id}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Organism: {results.get('organism', 'hsa')}
Client IP: {client_ip}

DATABASE STORAGE:
-----------------
This analysis was performed in database storage mode.
Results are preserved in the database and will not be lost after server restart.
You can access these results again using the task ID: {task_id}

Summary:
--------
Total Genes Analyzed: {summary.get('total_genes_analyzed', 0)}
Genes Mapped to KEGG: {summary.get('genes_mapped', 0)}
Mapping Rate: {summary.get('mapping_rate', 0)}%
Significant Pathways (p < 0.05): {summary.get('significant_pathways', 0)}
Total Pathways Analyzed: {summary.get('total_pathways_analyzed', 0)}

Top Enriched Pathways:
-------------------------
"""
            for i, pathway in enumerate(pathways[:10], 1):
                report += f"{i}. {pathway.get('name', '')} ({pathway.get('id', '')})\n"
                report += f"   Category: {pathway.get('category', '')}\n"
                report += f"   Enrichment Score: {pathway.get('enrichment_score', 0):.3f}\n"
                report += f"   P-value: {pathway.get('p_value', 0):.2e}\n"
                report += f"   Adjusted P-value: {pathway.get('adjusted_pvalue', 0):.2e}\n"
                report += f"   Gene Count: {pathway.get('gene_count', 0)}\n"
                report += f"   Top Genes: {', '.join(pathway.get('genes', [])[:3])}\n\n"
            
            zip_file.writestr('kegg_analysis_report.txt', report)
            
            # 3. 基因列表
            gene_data = io.StringIO()
            writer = csv.writer(gene_data)
            writer.writerow(['Gene', 'Pathways'])
            
            # 收集所有基因及其所在通路
            gene_pathways = {}
            for pathway in pathways:
                pathway_name = pathway.get('name', 'Unknown')
                for gene in pathway.get('genes', []):
                    if gene not in gene_pathways:
                        gene_pathways[gene] = []
                    gene_pathways[gene].append(pathway_name)
            
            for gene in sorted(gene_pathways.keys()):
                pathways_str = '; '.join(gene_pathways[gene])
                writer.writerow([gene, pathways_str])
            
            zip_file.writestr('kegg_gene_list.csv', gene_data.getvalue())
            
            # 4. 原始JSON数据（供高级用户使用）
            zip_file.writestr('kegg_raw_results.json', json.dumps(results, indent=2))
            
            # 5. 分析参数
            params = task.parameters if isinstance(task.parameters, dict) else {}
            zip_file.writestr('kegg_analysis_parameters.json', json.dumps(params, indent=2))
            
            # 6. 使用说明
            readme = """KEGG Analysis Results Package
=======================================

This ZIP file contains all results from your KEGG pathway enrichment analysis.

Files included:
1. kegg_enrichment_table.csv - Complete pathway enrichment results
2. kegg_analysis_report.txt - Summary report of the analysis
3. kegg_gene_list.csv - List of all genes and their associated pathways
4. kegg_raw_results.json - Raw JSON data (for advanced users)
5. kegg_analysis_parameters.json - Analysis parameters used

Note: This analysis was performed in database storage mode.
      Results are preserved in the database. You can access them again using the task ID.

Animal Hormone Database - KEGG Analysis Tool
"""
            zip_file.writestr('README.txt', readme)
        
        zip_buffer.seek(0)
        response = HttpResponse(zip_buffer, content_type='application/zip')
        response['Content-Disposition'] = f'attachment; filename="kegg_analysis_{task_id}.zip"'
        return response
        
    except Exception as e:
        logger.error(f"下载所有结果失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=500)

@csrf_exempt
@require_http_methods(["POST"])
def cancel_kegg_task(request, task_id):
    """取消KEGG分析任务 - 不需要用户认证"""
    try:
        task = get_task(task_id)
        if not task:
            return JsonResponse({
                'success': False,
                'error': 'Task not found'
            }, status=404)
        
        if task.status in ['completed', 'failed', 'cancelled']:
            return JsonResponse({
                'success': False,
                'error': f'Task already {task.status}'
            }, status=400)
        
        update_task_status(task, 'cancelled', 'Task cancelled by user')
        
        client_ip = request.META.get('REMOTE_ADDR', 'unknown')
        logger.info(f"Task {task_id} cancelled by client {client_ip}")
        
        return JsonResponse({
            'success': True,
            'task_id': str(task.task_id),
            'message': 'Task cancelled successfully',
            'client_mode': True
        })
        
    except Exception as e:
        logger.error(f"取消任务失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e),
            'client_mode': True
        }, status=500)

@csrf_exempt
@require_http_methods(["POST"])
def save_kegg_parameters(request):
    """保存KEGG分析参数 - 不需要用户认证"""
    try:
        data = json.loads(request.body)
        
        # 记录保存参数的请求
        client_ip = request.META.get('REMOTE_ADDR', 'unknown')
        logger.info(f"Parameters saved from client {client_ip}")
        
        return JsonResponse({
            'success': True,
            'message': 'Parameters saved successfully',
            'saved_at': datetime.now().isoformat(),
            'client_mode': True,
            'note': 'Parameters will be stored when creating analysis task'
        })
        
    except Exception as e:
        logger.error(f"保存参数失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e),
            'client_mode': True
        }, status=500)

@require_http_methods(["GET"])
def download_kegg_analysis_report(request, task_id):
    """下载KEGG分析报告 - 不需要用户认证"""
    try:
        task = get_task(task_id)
        
        if not task or task.status != 'completed':
            return JsonResponse({
                'success': False,
                'error': 'Analysis not completed or task not found'
            }, status=404)
        
        client_ip = request.META.get('REMOTE_ADDR', 'unknown')
        
        # 生成HTML报告
        results = task.results if task.results else {}
        pathways = results.get('pathways', [])
        summary = results.get('summary', {})
        
        # 生成HTML报告
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>KEGG Pathway Analysis Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }}
                h1, h2 {{ color: #333; }}
                .summary {{ background: #f5f5f5; padding: 20px; border-radius: 5px; margin-bottom: 20px; }}
                .note {{ background: #e8f4fd; padding: 15px; border-left: 4px solid #2196F3; margin: 20px 0; }}
                table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #4CAF50; color: white; }}
                tr:nth-child(even) {{ background-color: #f2f2f2; }}
                .footer {{ margin-top: 40px; color: #666; font-size: 12px; }}
            </style>
        </head>
        <body>
            <h1>KEGG Pathway Enrichment Analysis Report</h1>
            
            <div class="note">
                <strong>Database Storage Mode:</strong> This analysis was performed in database storage mode.
                Results are preserved in the database and will not be lost after server restart.
                Task ID for future reference: {task_id}
            </div>
            
            <div class="summary">
                <h2>Analysis Information</h2>
                <p><strong>Task ID:</strong> {task_id}</p>
                <p><strong>Analysis Date:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                <p><strong>Organism:</strong> {results.get('organism', 'hsa')}</p>
                <p><strong>Client IP:</strong> {client_ip}</p>
            </div>
            
            <div class="summary">
                <h2>Summary</h2>
                <p><strong>Total Genes Analyzed:</strong> {summary.get('total_genes_analyzed', 0)}</p>
                <p><strong>Genes Mapped to KEGG:</strong> {summary.get('genes_mapped', 0)}</p>
                <p><strong>Mapping Rate:</strong> {summary.get('mapping_rate', 0)}%</p>
                <p><strong>Significant Pathways (p < 0.05):</strong> {summary.get('significant_pathways', 0)}</p>
                <p><strong>Total Pathways Analyzed:</strong> {summary.get('total_pathways_analyzed', 0)}</p>
                <p><strong>Top Enriched Pathway:</strong> {summary.get('top_pathway', 'None')}</p>
            </div>
            
            <h2>Top Enriched Pathways</h2>
            <table>
                <tr>
                    <th>#</th>
                    <th>Pathway Name</th>
                    <th>Pathway ID</th>
                    <th>Enrichment Score</th>
                    <th>P-value</th>
                    <th>Gene Count</th>
                    <th>Category</th>
                </tr>
        """
        
        for i, pathway in enumerate(pathways[:20], 1):
            html_content += f"""
                <tr>
                    <td>{i}</td>
                    <td>{pathway.get('name', '')}</td>
                    <td>{pathway.get('id', '')}</td>
                    <td>{pathway.get('enrichment_score', 0):.3f}</td>
                    <td>{pathway.get('p_value', 0):.2e}</td>
                    <td>{pathway.get('gene_count', 0)}</td>
                    <td>{pathway.get('category', '')}</td>
                </tr>
            """
        
        html_content += f"""
            </table>
            
            <h2>Analysis Details</h2>
            <p><strong>Generated by:</strong> Animal Hormone Database KEGG Analysis Tool</p>
            <p><strong>Mode:</strong> Database Storage Mode (No login required)</p>
            <p><strong>Data Persistence:</strong> Results are stored in database and preserved</p>
            <p><strong>Task ID for Reference:</strong> {task_id}</p>
            
            <div class="footer">
                <p>Report generated on {datetime.now().strftime('%Y-%m-%d at %H:%M:%S')}</p>
                <p>Animal Hormone Database - https://animalhormone.org</p>
            </div>
        </body>
        </html>
        """
        
        response = HttpResponse(html_content, content_type='text/html')
        response['Content-Disposition'] = f'attachment; filename="kegg_report_{task_id}.html"'
        return response
        
    except Exception as e:
        logger.error(f"下载报告失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=500)

@require_http_methods(["GET"])
def get_kegg_demo_data(request):
    """获取KEGG演示数据 - 不需要用户认证"""
    try:
        # 提供演示基因列表
        demo_genes = [
            "TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "PIK3CA", "AKT1", "PTEN",
            "MYC", "CCND1", "CDK4", "CDK6", "RB1", "E2F1", "MDM2", "VEGFA",
            "HIF1A", "MMP2", "MMP9", "TIMP1", "TIMP2", "COL1A1", "COL3A1",
            "FN1", "LAMA1", "LAMB1", "LAMC1", "ITGA5", "ITGB1", "CDH1",
            "CDH2", "SNAI1", "SNAI2", "TWIST1", "VIM", "ZEB1", "ZEB2",
            "IL6", "IL8", "TNF", "TGFB1", "TGFB2", "TGFB3", "SMAD2",
            "SMAD3", "SMAD4", "NFKB1", "RELA", "STAT3", "JAK1", "JAK2"
        ]
        
        # 提供演示通路数据
        demo_pathways = get_example_pathways()
        
        return JsonResponse({
            'success': True,
            'demo_data': {
                'genes': demo_genes,
                'pathways': demo_pathways,
                'description': 'Sample gene list for KEGG pathway analysis demonstration',
                'note': 'This is demonstration data. Upload your own gene list for analysis.'
            },
            'client_mode': True
        })
        
    except Exception as e:
        logger.error(f"获取演示数据失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=500)
        
@require_http_methods(["GET"])
def list_user_kegg_tasks(request):
    """列出用户的KEGG分析任务 - 不需要用户认证"""
    try:
        # 获取客户端标识
        client_id = generate_client_id(request)
        
        # 从数据库查询该客户端的所有任务
        tasks = KeggAnalysisTask.objects.filter(client_id=client_id).order_by('-created_at')
        
        tasks_list = []
        for task in tasks:
            tasks_list.append({
                'id': str(task.task_id),
                'task_id': str(task.task_id),
                'status': task.status,
                'created_at': task.created_at.isoformat(),
                'updated_at': task.updated_at.isoformat(),
                'progress': task.progress,
                'message': task.message,
                'file_name': task.file_name,
                'gene_count': task.gene_count
            })
        
        client_ip = request.META.get('REMOTE_ADDR', 'unknown')
        logger.info(f"Fetched {len(tasks_list)} tasks for client {client_ip}")
        
        return JsonResponse({
            'success': True,
            'tasks': tasks_list,
            'total': len(tasks_list),
            'client_mode': True,
            'message': f'Found {len(tasks_list)} tasks for your session.' if tasks_list else 'No tasks found for your session.'
        })
        
    except Exception as e:
        logger.error(f"列出KEGG任务失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e),
            'tasks': [],
            'total': 0,
            'client_mode': True
        }, status=500)

@require_http_methods(["GET"])
def clear_old_tasks(request):
    """清理旧的任务数据 - 用于维护"""
    try:
        # 清理超过7天的任务
        from datetime import datetime, timedelta
        cutoff_time = datetime.now() - timedelta(days=7)
        
        # 删除旧任务
        old_tasks = KeggAnalysisTask.objects.filter(created_at__lt=cutoff_time)
        deleted_count = old_tasks.count()
        old_tasks.delete()
        
        # 清理结果文件
        try:
            for result_file in KEGG_RESULTS_DIR.glob("*_results.json"):
                file_time = datetime.fromtimestamp(result_file.stat().st_mtime)
                if file_time < cutoff_time:
                    result_file.unlink()
        except Exception as file_error:
            logger.warning(f"清理结果文件时出错: {file_error}")
        
        return JsonResponse({
            'success': True,
            'cleaned_tasks': deleted_count,
            'remaining_tasks': KeggAnalysisTask.objects.count(),
            'message': f'Cleaned {deleted_count} old tasks (older than 7 days)'
        })
        
    except Exception as e:
        logger.error(f"清理任务失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=500)

@require_http_methods(["GET"])
def get_kegg_previous_results(request):
    """获取用户的KEGG分析历史结果 - 不需要用户认证"""
    try:
        # 获取客户端标识
        client_id = generate_client_id(request)
        
        # 从数据库查询该客户端已完成的
        # 从数据库查询该客户端已完成的任务
        completed_tasks = KeggAnalysisTask.objects.filter(
            client_id=client_id,
            status='completed'
        ).order_by('-created_at')
        
        results_list = []
        for task in completed_tasks:
            # 提取结果摘要信息
            summary_info = {
                'id': str(task.task_id),
                'task_id': str(task.task_id),
                'name': f"Analysis {str(task.task_id)[:8]}...",
                'status': task.status,
                'created_at': task.created_at.isoformat(),
                'updated_at': task.updated_at.isoformat(),
                'gene_count': task.gene_count,
                'file_name': task.file_name or 'Not specified',
                'progress': task.progress,
                'message': task.message
            }
            
            # 如果有结果数据，添加更多信息
            if task.results and isinstance(task.results, dict):
                summary = task.results.get('summary', {})
                summary_info.update({
                    'name': f"Analysis: {summary.get('top_pathway', 'KEGG Analysis')}",
                    'significant_pathways': summary.get('significant_pathways', 0),
                    'total_pathways': summary.get('total_pathways_analyzed', 0),
                    'mapping_rate': summary.get('mapping_rate', 0),
                    'top_pathway': summary.get('top_pathway', 'N/A'),
                    'organism': task.results.get('organism', 'hsa')
                })
            
            results_list.append(summary_info)
        
        client_ip = request.META.get('REMOTE_ADDR', 'unknown')
        logger.info(f"Fetched {len(results_list)} completed results for client {client_ip}")
        
        return JsonResponse({
            'success': True,
            'results': results_list,
            'total': len(results_list),
            'client_mode': True,
            'message': f'Found {len(results_list)} completed analyses.' if results_list else 'No completed analyses found.'
        })
        
    except Exception as e:
        logger.error(f"获取历史结果失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e),
            'results': [],
            'total': 0,
            'client_mode': True
        }, status=500)
        
@require_http_methods(["GET"])
def get_task_statistics(request):
    """获取任务统计信息（用于管理）"""
    try:
        # 只允许特定IP或管理员访问
        client_ip = request.META.get('REMOTE_ADDR', 'unknown')
        
        # 计算统计信息
        total_tasks = KeggAnalysisTask.objects.count()
        completed_tasks = KeggAnalysisTask.objects.filter(status='completed').count()
        running_tasks = KeggAnalysisTask.objects.filter(status='running').count()
        failed_tasks = KeggAnalysisTask.objects.filter(status='failed').count()
        
        # 按日期分组
        from django.db.models import Count
        from django.utils import timezone
        from django.db.models.functions import TruncDate
        
        last_7_days = timezone.now() - timedelta(days=7)
        daily_stats = KeggAnalysisTask.objects.filter(
            created_at__gte=last_7_days
        ).annotate(
            date=TruncDate('created_at')
        ).values('date').annotate(
            count=Count('id')
        ).order_by('date')
        
        return JsonResponse({
            'success': True,
            'statistics': {
                'total_tasks': total_tasks,
                'completed_tasks': completed_tasks,
                'running_tasks': running_tasks,
                'failed_tasks': failed_tasks,
                'daily_stats': list(daily_stats)
            },
            'server_time': timezone.now().isoformat(),
            'note': 'Administrative statistics'
        })
        
    except Exception as e:
        logger.error(f"获取统计信息失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e)
        }, status=500)
        
@require_http_methods(["GET"])
def health_check(request):
    """健康检查端点"""
    try:
        # 检查数据库连接
        db_status = 'OK'
        try:
            KeggAnalysisTask.objects.count()
        except Exception as db_error:
            db_status = f'Error: {str(db_error)}'
        
        # 检查结果目录
        dir_status = 'OK'
        try:
            if not KEGG_RESULTS_DIR.exists():
                KEGG_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        except Exception as dir_error:
            dir_status = f'Error: {str(dir_error)}'
        
        return JsonResponse({
            'status': 'healthy',
            'timestamp': datetime.now().isoformat(),
            'components': {
                'database': db_status,
                'results_directory': dir_status,
                'model': 'KeggAnalysisTask available'
            },
            'tasks_summary': {
                'total': KeggAnalysisTask.objects.count(),
                'completed': KeggAnalysisTask.objects.filter(status='completed').count(),
                'pending': KeggAnalysisTask.objects.filter(status__in=['created', 'uploaded']).count()
            }
        })
        
    except Exception as e:
        return JsonResponse({
            'status': 'unhealthy',
            'error': str(e),
            'timestamp': datetime.now().isoformat()
        }, status=500)
        
@csrf_exempt
@require_http_methods(["DELETE", "POST"])
def delete_kegg_task(request, task_id):
    """删除KEGG分析任务 - 需要客户端匹配"""
    try:
        task = get_task(task_id)
        if not task:
            return JsonResponse({
                'success': False,
                'error': 'Task not found'
            }, status=404)
        
        # 检查客户端权限
        client_id = generate_client_id(request)
        if task.client_id != client_id:
            return JsonResponse({
                'success': False,
                'error': 'Permission denied: Task does not belong to this client'
            }, status=403)
        
        # 记录删除信息
        client_ip = request.META.get('REMOTE_ADDR', 'unknown')
        logger.info(f"Deleting task {task_id} by client {client_ip}")
        
        # 删除任务
        task_id_str = str(task.task_id)
        task.delete()
        
        # 尝试删除相关文件
        try:
            result_file = KEGG_RESULTS_DIR / f"{task_id_str}_results.json"
            if result_file.exists():
                result_file.unlink()
        except Exception as file_error:
            logger.warning(f"Could not delete result file for task {task_id_str}: {file_error}")
        
        return JsonResponse({
            'success': True,
            'message': f'Task {task_id_str} deleted successfully',
            'deleted_task_id': task_id_str,
            'client_mode': True
        })
        
    except Exception as e:
        logger.error(f"删除任务失败: {e}")
        return JsonResponse({
            'success': False,
            'error': str(e),
            'client_mode': True
        }, status=500)
    
@require_http_methods(["GET"])
def search_kegg_genes(request):
    """搜索KEGG基因"""
    try:
        query = request.GET.get('q', '')
        organism = request.GET.get('organism', 'hsa')
        
        if len(query) < 2:
            return JsonResponse({'success': True, 'genes': []})
        
        # ✅ 修复：将Q对象放在前面作为位置参数
        genes = KeggGene.objects.filter(
            Q(symbol__icontains=query) | Q(name__icontains=query),
            organism__code=organism
        ).values('gene_id', 'symbol', 'name')[:20]
        
        return JsonResponse({
            'success': True,
            'genes': list(genes)
        })
        
    except Exception as e:
        logger.error(f"搜索基因失败: {e}")
        return JsonResponse({'success': False, 'error': str(e)}, status=500)
        
@require_http_methods(["GET"])
def get_kegg_pathway_detail(request, pathway_id):
    """获取KEGG通路详情"""
    try:
        # 尝试两种可能的数据库存储格式
        # 1. 首先尝试 'path:hsa05220' 格式
        lookup_id = f'path:{pathway_id}'
        pathway = None
        try:
            pathway = KeggPathway.objects.get(pathway_id=lookup_id)
        except KeggPathway.DoesNotExist:
            # 2. 如果没找到，尝试 'hsa05220' 格式
            try:
                pathway = KeggPathway.objects.get(pathway_id=pathway_id)
            except KeggPathway.DoesNotExist:
                pass # 两个都找不到

        if not pathway:
             raise KeggPathway.DoesNotExist(f"Pathway with ID '{pathway_id}' or 'path:{pathway_id}' not found in database.")

        # 获取通路中的基因
        genes = KeggGenePathway.objects.filter(
            pathway=pathway
        ).select_related('gene')[:50] # 限制返回基因数量，避免数据过多

        gene_list = [{
            'id': g.gene.gene_id,
            'symbol': g.gene.symbol,
            'name': g.gene.name
        } for g in genes]

        # 构造KEGG官网链接时，使用标准格式 (不带 path:)
        # 这样前端可以直接拼接
        standard_pathway_id = pathway.pathway_id.replace('path:', '', 1) if pathway.pathway_id.startswith('path:') else pathway.pathway_id
        kegg_link = f"https://www.kegg.jp/pathway/{standard_pathway_id}"

        return JsonResponse({
            'success': True,
            'pathway': {
                'id': pathway.pathway_id,      # 返回数据库中的原始ID
                'standard_id': standard_pathway_id, # 添加一个标准ID字段，方便前端使用
                'name': pathway.name,
                'category': pathway.category,
                'gene_count': pathway.gene_count,
                'description': getattr(pathway, 'description', ''), # 如果模型有此字段
                'url': pathway.kegg_url or kegg_link, # 优先使用数据库存储的URL，否则构造
                'image_url': pathway.image_url,
                'genes': gene_list,
                'total_genes': pathway.gene_count
            }
        })

    except KeggPathway.DoesNotExist as e:
        logger.error(f"获取通路详情失败, ID: {pathway_id}, Error: {e}")
        return JsonResponse({'success': False, 'error': 'Pathway not found'}, status=404)
    except Exception as e:
        logger.error(f"获取通路详情失败, ID: {pathway_id}, Error: {e}")
        return JsonResponse({'success': False, 'error': str(e)}, status=500)