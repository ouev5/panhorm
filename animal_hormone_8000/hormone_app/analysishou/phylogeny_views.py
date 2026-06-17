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
from io import BytesIO, StringIO
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo, SeqIO, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# 导入您的现有模型
from ..models import (
    PhylogeneticTree,
    SequenceConservation,
    HormoneReceptorFull,
    HormoneReceptorInfo,
    HormoneRawData,
    HormoneRelatedGene
)

logger = logging.getLogger(__name__)

# 配置路径
PHYLOGENY_RESULTS_PATH = os.path.join(settings.MEDIA_ROOT, 'phylogeny_results')
os.makedirs(PHYLOGENY_RESULTS_PATH, exist_ok=True)

@csrf_exempt
@require_http_methods(["GET", "POST"])
def phylogeny_dashboard(request):
    """系统发育分析主页面"""
    return render(request, 'phylogeny_dashboard.html')

@csrf_exempt
@require_http_methods(["GET"])
def phylogeny_list_trees(request):
    """获取所有系统发育树列表"""
    try:
        trees = PhylogeneticTree.objects.all().order_by('-created_at')
        
        result = []
        for tree in trees:
            result.append({
                'id': tree.id,
                'hormone_family': tree.hormone_family,
                'tree_newick': tree.tree_newick[:200] + '...' if len(tree.tree_newick) > 200 else tree.tree_newick,
                'species_count': tree.species_count,
                'created_at': tree.created_at.isoformat() if tree.created_at else None
            })
        
        return JsonResponse({
            'success': True,
            'trees': result
        })
        
    except Exception as e:
        logger.error(f"Error listing trees: {str(e)}")
        return JsonResponse({'error': str(e)}, status=500)

@csrf_exempt
@require_http_methods(["GET"])
def phylogeny_get_tree(request, tree_id):
    """获取单个系统发育树详情"""
    try:
        tree = get_object_or_404(PhylogeneticTree, id=tree_id)
        
        return JsonResponse({
            'success': True,
            'tree': {
                'id': tree.id,
                'hormone_family': tree.hormone_family,
                'tree_newick': tree.tree_newick,
                'species_count': tree.species_count,
                'created_at': tree.created_at.isoformat() if tree.created_at else None
            }
        })
        
    except Exception as e:
        logger.error(f"Error getting tree: {str(e)}")
        return JsonResponse({'error': str(e)}, status=500)

@csrf_exempt
@require_http_methods(["GET"])
def phylogeny_get_conservation(request):
    """获取序列保守性分析结果"""
    try:
        hormone_id = request.GET.get('hormone_id')
        species = request.GET.get('species')
        
        queryset = SequenceConservation.objects.select_related('hormone').all()
        
        if hormone_id:
            queryset = queryset.filter(hormone_id=hormone_id)
        if species:
            queryset = queryset.filter(species__icontains=species)
        
        results = []
        for cons in queryset.order_by('-similarity_score')[:100]:
            hormone_name = ''
            if hasattr(cons.hormone, 'hormone_name'):
                hormone_name = cons.hormone.hormone_name
            elif hasattr(cons.hormone, 'receptor_name'):
                hormone_name = cons.hormone.receptor_name
            else:
                hormone_name = f"Hormone {cons.hormone.id}"
            
            results.append({
                'id': cons.id,
                'hormone': {
                    'id': cons.hormone.id,
                    'name': hormone_name
                },
                'species': cons.species,
                'similarity_score': cons.similarity_score,
                'alignment_length': cons.alignment_length
            })
        
        return JsonResponse({
            'success': True,
            'conservation': results
        })
        
    except Exception as e:
        logger.error(f"Error getting conservation: {str(e)}")
        return JsonResponse({'error': str(e)}, status=500)

@csrf_exempt
@require_http_methods(["GET"])
def phylogeny_get_families(request):
    """获取所有激素家族"""
    try:
        # 从 HormoneReceptorFull 获取激素家族
        families = set()
        family_counts = {}
        
        # 尝试从 HormoneReceptorFull 获取
        hormones = HormoneReceptorFull.objects.all()
        for h in hormones:
            # 从激素名称推断家族（例如："Insulin" -> "Insulin Family"）
            if hasattr(h, 'hormone_name') and h.hormone_name:
                name = h.hormone_name
                # 简单的家族推断逻辑
                if 'insulin' in name.lower():
                    family = 'Insulin Family'
                elif 'growth' in name.lower():
                    family = 'Growth Factor Family'
                elif 'steroid' in name.lower():
                    family = 'Steroid Hormone Family'
                elif 'thyroid' in name.lower():
                    family = 'Thyroid Hormone Family'
                else:
                    # 使用名称的第一个词作为家族
                    family = name.split()[0] + ' Family' if name.split() else 'Other'
                
                families.add(family)
                family_counts[family] = family_counts.get(family, 0) + 1
        
        # 也检查 HormoneReceptorInfo
        hormones2 = HormoneReceptorInfo.objects.all()
        for h in hormones2:
            if hasattr(h, 'hormone_name') and h.hormone_name:
                name = h.hormone_name
                if 'insulin' in name.lower():
                    family = 'Insulin Family'
                elif 'growth' in name.lower():
                    family = 'Growth Factor Family'
                elif 'steroid' in name.lower():
                    family = 'Steroid Hormone Family'
                elif 'thyroid' in name.lower():
                    family = 'Thyroid Hormone Family'
                else:
                    family = name.split()[0] + ' Family' if name.split() else 'Other'
                
                families.add(family)
                family_counts[family] = family_counts.get(family, 0) + 1
        
        family_list = [{'name': family, 'count': count} for family, count in family_counts.items()]
        
        return JsonResponse({
            'success': True,
            'families': family_list
        })
        
    except Exception as e:
        logger.error(f"Error getting families: {str(e)}")
        return JsonResponse({'error': str(e)}, status=500)

@csrf_exempt
@require_http_methods(["POST"])
def phylogeny_build_tree(request):
    """构建新的系统发育树"""
    try:
        data = json.loads(request.body)
        family_name = data.get('family_name')
        
        if not family_name:
            return JsonResponse({'error': 'Family name is required'}, status=400)
        
        # 获取该家族的激素序列
        sequences = []
        
        # 从 HormoneReceptorFull 获取序列
        hormones = HormoneReceptorFull.objects.all()
        for h in hormones:
            if hasattr(h, 'hormone_name') and h.hormone_name:
                # 检查是否属于所选家族
                is_family = False
                if family_name == 'Insulin Family' and 'insulin' in h.hormone_name.lower():
                    is_family = True
                elif family_name == 'Growth Factor Family' and 'growth' in h.hormone_name.lower():
                    is_family = True
                elif family_name == 'Steroid Hormone Family' and 'steroid' in h.hormone_name.lower():
                    is_family = True
                elif family_name == 'Thyroid Hormone Family' and 'thyroid' in h.hormone_name.lower():
                    is_family = True
                elif h.hormone_name.split()[0] + ' Family' == family_name:
                    is_family = True
                
                if is_family and hasattr(h, 'hormone_coding_genes_sequence') and h.hormone_coding_genes_sequence:
                    sequences.append({
                        'id': h.id,
                        'name': h.hormone_name,
                        'species': getattr(h, 'hormone_species_name', 'Unknown'),
                        'sequence': h.hormone_coding_genes_sequence
                    })
        
        # 也从 HormoneReceptorInfo 获取
        hormones2 = HormoneReceptorInfo.objects.all()
        for h in hormones2:
            if hasattr(h, 'hormone_name') and h.hormone_name:
                is_family = False
                if family_name == 'Insulin Family' and 'insulin' in h.hormone_name.lower():
                    is_family = True
                elif family_name == 'Growth Factor Family' and 'growth' in h.hormone_name.lower():
                    is_family = True
                elif family_name == 'Steroid Hormone Family' and 'steroid' in h.hormone_name.lower():
                    is_family = True
                elif family_name == 'Thyroid Hormone Family' and 'thyroid' in h.hormone_name.lower():
                    is_family = True
                elif h.hormone_name.split()[0] + ' Family' == family_name:
                    is_family = True
                
                if is_family and hasattr(h, 'receptor_coding_genes_sequence') and h.receptor_coding_genes_sequence:
                    sequences.append({
                        'id': h.id,
                        'name': h.hormone_name,
                        'species': getattr(h, 'receptor_species_name', 'Unknown'),
                        'sequence': h.receptor_coding_genes_sequence
                    })
        
        if len(sequences) < 3:
            return JsonResponse({
                'error': f'Need at least 3 sequences, only {len(sequences)} found'
            }, status=400)
        
        # 创建临时FASTA文件
        fasta_path = os.path.join(PHYLOGENY_RESULTS_PATH, f'{uuid.uuid4()}.fasta')
        with open(fasta_path, 'w') as f:
            for seq in sequences:
                f.write(f'>{seq["id"]}|{seq["species"]}|{seq["name"]}\n')
                f.write(f'{seq["sequence"]}\n')
        
        # 读取序列
        records = list(SeqIO.parse(fasta_path, 'fasta'))
        
        # ========== 修复部分开始 ==========
        # 计算距离矩阵
        n = len(records)
        
        # 创建距离矩阵
        distance_matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i+1, n):
                seq1 = str(records[i].seq)
                seq2 = str(records[j].seq)
                
                # 简单的相似度计算
                min_len = min(len(seq1), len(seq2))
                if min_len == 0:
                    similarity = 0
                else:
                    matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a == b)
                    similarity = matches / min_len
                
                distance = 1 - similarity
                distance_matrix[i][j] = distance
                distance_matrix[j][i] = distance
        
        # 构建下三角矩阵（_DistanceMatrix 需要的格式）
        from Bio.Phylo.TreeConstruction import _DistanceMatrix
        
        # 关键修复：正确构建下三角矩阵
        lower_triangular = []
        for i in range(1, n):  # 从第1行开始（因为第0行没有元素）
            row = []
            for j in range(i):  # j < i
                row.append(float(distance_matrix[i][j]))
            lower_triangular.append(row)
        
        # 打印调试信息
        logger.info(f"下三角矩阵行数: {len(lower_triangular)}")
        for idx, row in enumerate(lower_triangular):
            logger.info(f"行 {idx+1}: {row}")
        
        names = [str(r.id) for r in records]
        
        # 创建距离矩阵对象 - 注意：_DistanceMatrix 期望的是所有行（包括空的第一行）
        # 实际上它内部会处理，我们传入的就是下三角矩阵
        dm = _DistanceMatrix(names=names, matrix=lower_triangular)
        # ========== 修复部分结束 ==========
        
        # 构建树
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        
        # 转换为Newick格式
        tree.ladderize()
        newick_io = StringIO()
        Phylo.write(tree, newick_io, 'newick')
        newick = newick_io.getvalue()
        
        # 保存到数据库
        phy_tree, created = PhylogeneticTree.objects.update_or_create(
            hormone_family=family_name,
            defaults={
                'tree_newick': newick,
                'species_count': len(sequences)
            }
        )
        
        # 计算保守性分数并保存
        for i, record in enumerate(records):
            try:
                hormone_id = int(record.id.split('|')[0])
                species_name = record.id.split('|')[1] if len(record.id.split('|')) > 1 else 'Unknown'
                
                # 计算与其他序列的平均相似度
                similarities = []
                for j, other in enumerate(records):
                    if i != j:
                        seq1 = str(record.seq)
                        seq2 = str(other.seq)
                        min_len = min(len(seq1), len(seq2))
                        if min_len > 0:
                            matches = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a == b)
                            similarity = matches / min_len
                            similarities.append(similarity)
                
                avg_similarity = sum(similarities) / len(similarities) if similarities else 0
                
                # 查找对应的激素模型
                hormone_model = None
                try:
                    hormone_model = HormoneReceptorFull.objects.get(id=hormone_id)
                except HormoneReceptorFull.DoesNotExist:
                    try:
                        hormone_model = HormoneReceptorInfo.objects.get(id=hormone_id)
                    except HormoneReceptorInfo.DoesNotExist:
                        pass
                
                if hormone_model:
                    SequenceConservation.objects.update_or_create(
                        hormone=hormone_model,
                        species=species_name,
                        defaults={
                            'similarity_score': avg_similarity,
                            'alignment_length': len(record.seq)
                        }
                    )
            except (ValueError, IndexError) as e:
                logger.warning(f"Could not parse record ID: {record.id}")
                continue
        
        # 清理临时文件
        if os.path.exists(fasta_path):
            os.remove(fasta_path)
        
        return JsonResponse({
            'success': True,
            'tree_id': phy_tree.id,
            'message': f'Tree built successfully for {family_name} with {len(sequences)} sequences'
        })
        
    except Exception as e:
        logger.error(f"Error building tree: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return JsonResponse({'error': str(e)}, status=500)

@csrf_exempt
@require_http_methods(["GET"])
def phylogeny_download_tree(request, tree_id):
    """下载Newick格式的系统发育树"""
    try:
        tree = get_object_or_404(PhylogeneticTree, id=tree_id)
        
        response = FileResponse(
            BytesIO(tree.tree_newick.encode('utf-8')),
            content_type='text/plain'
        )
        filename = f"{tree.hormone_family.replace(' ', '_')}_tree_{timezone.now().strftime('%Y%m%d')}.nwk"
        response['Content-Disposition'] = f'attachment; filename="{filename}"'
        return response
        
    except Exception as e:
        logger.error(f"Error downloading tree: {str(e)}")
        raise Http404("Could not download tree file")

@csrf_exempt
@require_http_methods(["GET"])
def phylogeny_download_conservation_csv(request):
    """下载保守性分析结果CSV"""
    try:
        hormone_id = request.GET.get('hormone_id')
        
        queryset = SequenceConservation.objects.select_related('hormone').all()
        if hormone_id:
            queryset = queryset.filter(hormone_id=hormone_id)
        
        data = []
        for cons in queryset:
            hormone_name = ''
            if hasattr(cons.hormone, 'hormone_name'):
                hormone_name = cons.hormone.hormone_name
            elif hasattr(cons.hormone, 'receptor_name'):
                hormone_name = cons.hormone.receptor_name
            
            data.append({
                'Hormone': hormone_name,
                'Hormone ID': cons.hormone.id,
                'Species': cons.species,
                'Similarity Score': cons.similarity_score,
                'Alignment Length': cons.alignment_length
            })
        
        df = pd.DataFrame(data)
        
        csv_buffer = BytesIO()
        df.to_csv(csv_buffer, index=False, encoding='utf-8-sig')
        csv_buffer.seek(0)
        
        response = FileResponse(
            csv_buffer,
            content_type='text/csv'
        )
        filename = f"conservation_{timezone.now().strftime('%Y%m%d_%H%M%S')}.csv"
        response['Content-Disposition'] = f'attachment; filename="{filename}"'
        return response
        
    except Exception as e:
        logger.error(f"Error downloading conservation CSV: {str(e)}")
        raise Http404("Could not generate CSV file")

@csrf_exempt
@require_http_methods(["GET"])
def phylogeny_download_heatmap(request):
    """下载保守性热力图"""
    try:
        hormone_id = request.GET.get('hormone_id')
        
        queryset = SequenceConservation.objects.select_related('hormone').all()
        if hormone_id:
            queryset = queryset.filter(hormone_id=hormone_id)
        
        if queryset.count() == 0:
            raise Http404("No conservation data available")
        
        # 准备数据用于热力图
        data = []
        hormones = set()
        species = set()
        
        for cons in queryset:
            hormone_name = ''
            if hasattr(cons.hormone, 'hormone_name'):
                hormone_name = cons.hormone.hormone_name
            elif hasattr(cons.hormone, 'receptor_name'):
                hormone_name = cons.hormone.receptor_name
            
            hormones.add(hormone_name)
            species.add(cons.species)
            data.append({
                'hormone': hormone_name,
                'species': cons.species,
                'score': cons.similarity_score
            })
        
        # 创建矩阵
        hormone_list = sorted(list(hormones))
        species_list = sorted(list(species))
        
        matrix = np.zeros((len(hormone_list), len(species_list)))
        for item in data:
            i = hormone_list.index(item['hormone'])
            j = species_list.index(item['species'])
            matrix[i, j] = item['score']
        
        # 绘制热力图
        plt.figure(figsize=(max(10, len(species_list) * 0.5), max(8, len(hormone_list) * 0.4)))
        
        sns.heatmap(matrix, 
                   xticklabels=species_list, 
                   yticklabels=hormone_list,
                   cmap='YlOrRd',
                   annot=True,
                   fmt='.2f',
                   annot_kws={'size': 8},
                   cbar_kws={'label': 'Similarity Score'})
        
        plt.title('Sequence Conservation Heatmap Across Species', fontsize=14, fontweight='bold')
        plt.xlabel('Species', fontsize=12)
        plt.ylabel('Hormone', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        buf = BytesIO()
        plt.savefig(buf, format='png', dpi=150, bbox_inches='tight')
        plt.close()
        buf.seek(0)
        
        response = FileResponse(buf, content_type='image/png')
        filename = f"conservation_heatmap_{timezone.now().strftime('%Y%m%d_%H%M%S')}.png"
        response['Content-Disposition'] = f'attachment; filename="{filename}"'
        return response
        
    except Exception as e:
        logger.error(f"Error generating heatmap: {str(e)}")
        raise Http404("Could not generate heatmap")