import os
import json
import uuid
import logging
import pandas as pd
import numpy as np
from datetime import datetime
from django.shortcuts import render, get_object_or_404
from django.http import JsonResponse, FileResponse, Http404
from django.conf import settings
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_http_methods
from django.utils import timezone
from django.db import models
from django.db.models import Q
from io import BytesIO
import matplotlib.pyplot as plt
import seaborn as sns

# 导入您的现有模型
from ..models import (
    GNPrediction, 
    HormoneReceptorFull,  # 肽类激素受体信息
    HormoneReceptorInfo,  # 非肽类激素受体信息
)

logger = logging.getLogger(__name__)

# 配置路径
GNN_RESULTS_PATH = os.path.join(settings.MEDIA_ROOT, 'gnn_results')
os.makedirs(GNN_RESULTS_PATH, exist_ok=True)

@csrf_exempt
@require_http_methods(["GET", "POST"])
def gnn_dashboard(request):
    """GNN预测主页面"""
    return render(request, 'gnn_dashboard.html')

@csrf_exempt
@require_http_methods(["GET"])
def gnn_get_predictions(request):
    """获取GNN预测结果列表"""
    try:
        # 获取查询参数
        hormone_id = request.GET.get('hormone_id')
        receptor_id = request.GET.get('receptor_id')
        confidence = request.GET.get('confidence')
        hormone_type = request.GET.get('hormone_type')
        receptor_type = request.GET.get('receptor_type')
        limit = int(request.GET.get('limit', 100))
        
        # 构建查询
        queryset = GNPrediction.objects.all()
        
        # 应用过滤器
        if hormone_id:
            queryset = queryset.filter(
                Q(hormone_peptide_id=hormone_id) | 
                Q(hormone_non_peptide_id=hormone_id)
            )
        if receptor_id:
            queryset = queryset.filter(
                Q(receptor_peptide_id=receptor_id) | 
                Q(receptor_non_peptide_id=receptor_id)
            )
        if confidence:
            queryset = queryset.filter(confidence_level=confidence)
        if hormone_type:
            queryset = queryset.filter(hormone_type=hormone_type)
        if receptor_type:
            queryset = queryset.filter(receptor_type=receptor_type)
        
        # 按分数排序并限制数量
        queryset = queryset.order_by('-prediction_score')[:limit]
        
        # 序列化结果
        predictions = []
        for pred in queryset:
            # 获取激素信息
            hormone_info = get_hormone_info(pred)
            
            # 获取受体信息
            receptor_info = get_receptor_info(pred)
            
            predictions.append({
                'id': pred.id,
                'hormone': {
                    'id': hormone_info['id'],
                    'name': hormone_info['name'],
                    'type': pred.hormone_type,
                    'hormone_class': '肽类激素' if pred.hormone_type == 'peptide' else '非肽类激素',
                    'species': hormone_info['species']
                },
                'receptor': {
                    'id': receptor_info['id'],
                    'name': receptor_info['name'],
                    'type': pred.receptor_type,
                    'gene_symbol': receptor_info['gene_symbol'],
                    'species': receptor_info['species']
                },
                'prediction_score': pred.prediction_score,
                'confidence_level': pred.confidence_level,
                'model_version': pred.model_version,
                'created_at': pred.created_at.isoformat() if pred.created_at else None
            })
        
        # 获取统计信息
        stats = {
            'total_predictions': GNPrediction.objects.count(),
            'high_confidence': GNPrediction.objects.filter(confidence_level='high').count(),
            'medium_confidence': GNPrediction.objects.filter(confidence_level='medium').count(),
            'low_confidence': GNPrediction.objects.filter(confidence_level='low').count(),
            'avg_score': GNPrediction.objects.aggregate(avg=models.Avg('prediction_score'))['avg'] or 0,
            'peptide_count': GNPrediction.objects.filter(hormone_type='peptide').count(),
            'non_peptide_count': GNPrediction.objects.filter(hormone_type='non_peptide').count()
        }
        
        return JsonResponse({
            'success': True,
            'predictions': predictions,
            'stats': stats
        })
        
    except Exception as e:
        logger.error(f"Error getting GNN predictions: {str(e)}")
        return JsonResponse({'error': str(e)}, status=500)

def get_hormone_info(pred):
    """辅助函数：获取激素信息"""
    result = {
        'id': None,
        'name': '未知激素',
        'species': '',
        'uniprot': '',
        'sequence': '',
        'genes': []
    }
    
    if pred.hormone_type == 'peptide' and pred.hormone_peptide_id:
        try:
            hormone = HormoneReceptorFull.objects.get(id=pred.hormone_peptide_id)
            result['id'] = hormone.id
            result['name'] = hormone.hormone_name or f"肽类激素-{pred.hormone_peptide_id}"
            result['species'] = hormone.hormone_species_name or ''
            result['uniprot'] = hormone.hormone_uniprot_id or ''
            result['sequence'] = hormone.hormone_coding_genes_sequence or ''
            if hormone.hormone_coding_genes:
                result['genes'].append(hormone.hormone_coding_genes)
        except HormoneReceptorFull.DoesNotExist:
            result['name'] = f"肽类激素(已删除)"
            result['id'] = pred.hormone_peptide_id
            
    elif pred.hormone_type == 'non_peptide' and pred.hormone_non_peptide_id:
        try:
            hormone = HormoneReceptorInfo.objects.get(id=pred.hormone_non_peptide_id)
            result['id'] = hormone.id
            result['name'] = hormone.hormone_name or f"非肽类激素-{pred.hormone_non_peptide_id}"
            result['species'] = hormone.receptor_species_name or ''
            result['uniprot'] = hormone.receptor_uniprot_id or ''
            result['sequence'] = hormone.receptor_coding_genes_sequence or ''
            if hormone.receptor_coding_genes:
                result['genes'].append(hormone.receptor_coding_genes)
        except HormoneReceptorInfo.DoesNotExist:
            result['name'] = f"非肽类激素(已删除)"
            result['id'] = pred.hormone_non_peptide_id
    
    return result

def get_receptor_info(pred):
    """辅助函数：获取受体信息"""
    result = {
        'id': None,
        'name': '未知受体',
        'species': '',
        'uniprot': '',
        'sequence': '',
        'gene_symbol': '',
        'genes': []
    }
    
    if pred.receptor_type == 'peptide' and pred.receptor_peptide_id:
        try:
            receptor = HormoneReceptorFull.objects.get(id=pred.receptor_peptide_id)
            result['id'] = receptor.id
            result['name'] = receptor.receptor_name or f"肽类受体-{pred.receptor_peptide_id}"
            result['species'] = receptor.receptor_species_name or ''
            result['uniprot'] = receptor.receptor_uniprot_id or ''
            result['sequence'] = receptor.receptor_coding_genes_sequence or ''
            result['gene_symbol'] = receptor.receptor_coding_genes or ''
            if receptor.receptor_coding_genes:
                result['genes'].append(receptor.receptor_coding_genes)
        except HormoneReceptorFull.DoesNotExist:
            result['name'] = f"肽类受体(已删除)"
            result['id'] = pred.receptor_peptide_id
            
    elif pred.receptor_type == 'non_peptide' and pred.receptor_non_peptide_id:
        try:
            receptor = HormoneReceptorInfo.objects.get(id=pred.receptor_non_peptide_id)
            result['id'] = receptor.id
            result['name'] = receptor.receptor_name or f"非肽类受体-{pred.receptor_non_peptide_id}"
            result['species'] = receptor.receptor_species_name or ''
            result['uniprot'] = receptor.receptor_uniprot_id or ''
            result['sequence'] = receptor.receptor_coding_genes_sequence or ''
            result['gene_symbol'] = receptor.receptor_coding_genes or ''
            if receptor.receptor_coding_genes:
                result['genes'].append(receptor.receptor_coding_genes)
        except HormoneReceptorInfo.DoesNotExist:
            result['name'] = f"非肽类受体(已删除)"
            result['id'] = pred.receptor_non_peptide_id
    
    return result

@csrf_exempt
@require_http_methods(["GET"])
def gnn_get_prediction_detail(request, prediction_id):
    """获取单个预测详情"""
    try:
        prediction = get_object_or_404(GNPrediction, id=prediction_id)
        
        # 获取激素信息
        hormone_info = get_hormone_info(prediction)
        
        # 获取受体信息
        receptor_info = get_receptor_info(prediction)
        
        # 获取相关基因信息
        hormone_genes = hormone_info['genes']
        receptor_genes = receptor_info['genes']
        
        return JsonResponse({
            'success': True,
            'prediction': {
                'id': prediction.id,
                'hormone': {
                    'id': hormone_info['id'],
                    'name': hormone_info['name'],
                    'uniprot_id': hormone_info['uniprot'],
                    'type': prediction.hormone_type,
                    'species': hormone_info['species'],
                    'sequence': hormone_info['sequence'][:200] + '...' if hormone_info['sequence'] and len(hormone_info['sequence']) > 200 else hormone_info['sequence']
                },
                'receptor': {
                    'id': receptor_info['id'],
                    'name': receptor_info['name'],
                    'uniprot_id': receptor_info['uniprot'],
                    'gene_symbol': receptor_info['gene_symbol'],
                    'type': prediction.receptor_type,
                    'species': receptor_info['species'],
                    'sequence': receptor_info['sequence'][:200] + '...' if receptor_info['sequence'] and len(receptor_info['sequence']) > 200 else receptor_info['sequence']
                },
                'prediction_score': prediction.prediction_score,
                'confidence_level': prediction.confidence_level,
                'model_version': prediction.model_version,
                'created_at': prediction.created_at.isoformat() if prediction.created_at else None,
                'hormone_genes': [{'symbol': g, 'name': g} for g in hormone_genes],
                'receptor_genes': [{'symbol': g, 'name': g} for g in receptor_genes]
            }
        })
        
    except Exception as e:
        logger.error(f"Error getting prediction detail: {str(e)}")
        return JsonResponse({'error': str(e)}, status=500)

@csrf_exempt
@require_http_methods(["GET"])
def gnn_get_network_graph(request):
    """获取知识图谱网络数据"""
    try:
        hormone_id = request.GET.get('hormone_id')
        
        # 构建节点和边
        nodes = []
        edges = []
        node_ids = set()
        
        # 获取预测数据
        predictions = GNPrediction.objects.all()
        if hormone_id:
            predictions = predictions.filter(
                Q(hormone_peptide_id=hormone_id) | 
                Q(hormone_non_peptide_id=hormone_id)
            )
        
        predictions = predictions[:50]  # 限制数量避免图过大
        
        for pred in predictions:
            # 获取激素信息
            hormone_info = get_hormone_info(pred)
            
            # 添加激素节点
            hormone_key = f'h_{hormone_info["id"]}'
            if hormone_key not in node_ids and hormone_info['id']:
                hormone_color = '#165DFF' if pred.hormone_type == 'peptide' else '#36CFC9'
                
                nodes.append({
                    'id': hormone_key,
                    'label': hormone_info['name'][:20],
                    'type': 'hormone',
                    'group': pred.hormone_type,
                    'size': 20,
                    'color': hormone_color
                })
                node_ids.add(hormone_key)
            
            # 获取受体信息
            receptor_info = get_receptor_info(pred)
            
            # 添加受体节点
            receptor_key = f'r_{receptor_info["id"]}'
            if receptor_key not in node_ids and receptor_info['id']:
                receptor_color = '#52C41A' if pred.receptor_type == 'peptide' else '#FAAD14'
                
                nodes.append({
                    'id': receptor_key,
                    'label': receptor_info['name'][:20],
                    'type': 'receptor',
                    'group': pred.receptor_type,
                    'size': 15,
                    'color': receptor_color
                })
                node_ids.add(receptor_key)
            
            # 添加边
            if hormone_info['id'] and receptor_info['id']:
                edges.append({
                    'source': hormone_key,
                    'target': receptor_key,
                    'weight': pred.prediction_score,
                    'confidence': pred.confidence_level,
                    'type': 'predicted'
                })
        
        return JsonResponse({
            'success': True,
            'nodes': nodes,
            'edges': edges
        })
        
    except Exception as e:
        logger.error(f"Error generating network graph: {str(e)}")
        return JsonResponse({'error': str(e)}, status=500)

@csrf_exempt
@require_http_methods(["GET"])
def gnn_download_predictions_csv(request):
    """下载预测结果为CSV"""
    try:
        # 创建CSV数据
        predictions = GNPrediction.objects.all()
        
        data = []
        for pred in predictions:
            hormone_info = get_hormone_info(pred)
            receptor_info = get_receptor_info(pred)
            
            data.append({
                'Hormone': hormone_info['name'],
                'Hormone ID': hormone_info['id'],
                'Hormone Type': '肽类' if pred.hormone_type == 'peptide' else '非肽类',
                'Hormone Species': hormone_info['species'],
                'Receptor': receptor_info['name'],
                'Receptor ID': receptor_info['id'],
                'Receptor Type': '肽类' if pred.receptor_type == 'peptide' else '非肽类',
                'Receptor Species': receptor_info['species'],
                'Gene Symbol': receptor_info['gene_symbol'],
                'Prediction Score': pred.prediction_score,
                'Confidence Level': pred.confidence_level,
                'Model Version': pred.model_version,
                'Prediction Date': pred.created_at.strftime('%Y-%m-%d') if pred.created_at else ''
            })
        
        df = pd.DataFrame(data)
        
        # 创建CSV文件
        csv_buffer = BytesIO()
        df.to_csv(csv_buffer, index=False, encoding='utf-8-sig')
        csv_buffer.seek(0)
        
        response = FileResponse(
            csv_buffer,
            content_type='text/csv'
        )
        filename = f"gnn_predictions_{timezone.now().strftime('%Y%m%d_%H%M%S')}.csv"
        response['Content-Disposition'] = f'attachment; filename="{filename}"'
        return response
        
    except Exception as e:
        logger.error(f"Error downloading CSV: {str(e)}")
        raise Http404("Could not generate CSV file")

@csrf_exempt
@require_http_methods(["GET"])
def gnn_download_visualization(request, viz_type):
    """下载可视化图表"""
    try:
        predictions = GNPrediction.objects.all()
        
        if predictions.count() == 0:
            raise Http404("No predictions available for visualization")
        
        plt.figure(figsize=(12, 8))
        plt.style.use('seaborn-v0_8-darkgrid')
        
        if viz_type == 'score_distribution':
            # 分数分布直方图
            scores = [p.prediction_score for p in predictions]
            plt.hist(scores, bins=30, edgecolor='black', alpha=0.7, color='#722ED1')
            plt.title('Prediction Score Distribution', fontsize=14, fontweight='bold')
            plt.xlabel('Score', fontsize=12)
            plt.ylabel('Frequency', fontsize=12)
            plt.grid(True, alpha=0.3)
            
        elif viz_type == 'confidence_pie':
            # 置信度饼图
            confidence_counts = [
                predictions.filter(confidence_level='high').count(),
                predictions.filter(confidence_level='medium').count(),
                predictions.filter(confidence_level='low').count()
            ]
            labels = ['High', 'Medium', 'Low']
            colors = ['#52C41A', '#FAAD14', '#FF4D4F']
            
            # 过滤掉0值
            non_zero = [(count, label, color) for count, label, color in zip(confidence_counts, labels, colors) if count > 0]
            if non_zero:
                counts, plot_labels, plot_colors = zip(*non_zero)
                plt.pie(
                    counts,
                    labels=[f'{l} ({c})' for l, c in zip(plot_labels, counts)],
                    colors=plot_colors,
                    autopct='%1.1f%%',
                    startangle=90,
                    wedgeprops={'edgecolor': 'white', 'linewidth': 1}
                )
                plt.title('Confidence Level Distribution', fontsize=14, fontweight='bold')
                plt.axis('equal')
            
        elif viz_type == 'top_predictions':
            # 前20个预测的条形图
            top_preds = predictions.order_by('-prediction_score')[:20]
            names = []
            scores = []
            
            for p in top_preds:
                hormone_info = get_hormone_info(p)
                receptor_info = get_receptor_info(p)
                
                label = f"{hormone_info['name']} → {receptor_info['name']}"
                if len(label) > 40:
                    label = label[:37] + "..."
                names.append(label)
                scores.append(p.prediction_score)
            
            y_pos = range(len(names))
            plt.barh(y_pos, scores, color='#165DFF')
            plt.yticks(y_pos, names, fontsize=10)
            plt.xlabel('Prediction Score', fontsize=12)
            plt.title('Top 20 Predicted Interactions', fontsize=14, fontweight='bold')
            plt.gca().invert_yaxis()
            plt.tight_layout()
        
        # 保存到缓冲区
        buf = BytesIO()
        plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
        plt.close()
        buf.seek(0)
        
        response = FileResponse(buf, content_type='image/png')
        filename = f"gnn_{viz_type}_{timezone.now().strftime('%Y%m%d_%H%M%S')}.png"
        response['Content-Disposition'] = f'attachment; filename="{filename}"'
        return response
        
    except Exception as e:
        logger.error(f"Error generating visualization: {str(e)}")
        raise Http404("Could not generate visualization")