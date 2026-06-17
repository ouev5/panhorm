# services/kegg_analysis.py
import math
import numpy as np
from scipy import stats
from collections import defaultdict
from django.db.models import Count, Q
from ..models import KeggGene, KeggPathway, KeggGenePathway, KeggOrganism

class KeggAnalysisService:
    def __init__(self, organism='hsa'):
        self.organism = organism
        self.organism_obj = KeggOrganism.objects.get(code=organism)
    
    def perform_enrichment_analysis(self, gene_symbols, method='hypergeometric'):
        """执行KEGG通路富集分析"""
        
        # 1. 映射基因符号到KEGG基因
        mapped_genes = self.map_gene_symbols(gene_symbols)
        
        if not mapped_genes:
            return {'error': 'No genes mapped to KEGG database'}
        
        # 2. 获取背景信息
        total_genes_in_genome = KeggGene.objects.filter(organism=self.organism_obj).count()
        
        # 3. 获取所有通路
        all_pathways = KeggPathway.objects.filter(
            gene_relations__gene__organism=self.organism_obj
        ).distinct()
        
        # 4. 执行富集分析
        results = []
        for pathway in all_pathways:
            # 获取通路中的基因
            pathway_genes = set(
                KeggGene.objects.filter(
                    organism=self.organism_obj,
                    pathway_relations__pathway=pathway
                ).values_list('symbol', flat=True)
            )
            
            # 计算重叠
            overlap_genes = set(mapped_genes) & pathway_genes
            overlap_count = len(overlap_genes)
            
            if overlap_count == 0:
                continue
            
            # 计算统计量
            stats_result = self.calculate_enrichment_stats(
                total_genes_in_genome,
                len(pathway_genes),
                len(mapped_genes),
                overlap_count,
                method
            )
            
            # 构建结果
            result = {
                'pathway_id': pathway.pathway_id,
                'pathway_name': pathway.name,
                'category': pathway.category,
                'pathway_gene_count': len(pathway_genes),
                'input_gene_count': len(mapped_genes),
                'overlap_count': overlap_count,
                'overlap_genes': list(overlap_genes),
                'overlap_ratio': overlap_count / len(mapped_genes),
                'pathway_url': pathway.kegg_url or f"https://www.kegg.jp/pathway/{pathway.pathway_id}",
                'pathway_image': pathway.image_url or f"https://www.kegg.jp/kegg/pathway/{self.organism}/{pathway.pathway_id}.png",
            }
            result.update(stats_result)
            
            results.append(result)
        
        # 5. 多重检验校正
        if results:
            results = self.adjust_pvalues(results)
            results.sort(key=lambda x: x['p_value'])
        
        return {
            'organism': self.organism,
            'input_gene_count': len(gene_symbols),
            'mapped_gene_count': len(mapped_genes),
            'mapping_rate': len(mapped_genes) / len(gene_symbols) if gene_symbols else 0,
            'significant_pathways': len([r for r in results if r['adjusted_pvalue'] < 0.05]),
            'total_pathways_analyzed': len(results),
            'results': results[:100]  # 只返回前100个
        }
    
    def map_gene_symbols(self, gene_symbols):
        """映射基因符号到KEGG数据库"""
        # 清理基因符号
        cleaned_symbols = [g.upper().strip() for g in gene_symbols if g.strip()]
        
        # 查询数据库
        mapped_genes = KeggGene.objects.filter(
            organism=self.organism_obj,
            symbol__in=cleaned_symbols
        ).values_list('symbol', flat=True)
        
        # 如果没有找到，尝试别名映射
        if len(mapped_genes) < len(cleaned_symbols) * 0.5:  # 映射率低于50%
            # 这里可以添加别名映射逻辑
            pass
        
        return list(mapped_genes)
    
    def calculate_enrichment_stats(self, N, K, n, k, method='hypergeometric'):
        """计算富集统计量"""
        if method == 'hypergeometric':
            # 超几何检验
            p_value = stats.hypergeom.sf(k-1, N, K, n)
            
            # 富集分数
            if k == 0:
                enrichment_score = 0
            else:
                expected = n * K / N
                enrichment_score = k / expected if expected > 0 else 0
            
            return {
                'p_value': p_value,
                'enrichment_score': enrichment_score,
                'odds_ratio': (k / n) / ((K - k) / (N - n)) if (K - k) > 0 else float('inf'),
                'confidence_interval': self.calculate_confidence_interval(k, K, n, N)
            }
        
        return {'p_value': 1.0, 'enrichment_score': 0}
    
    def adjust_pvalues(self, results, method='fdr_bh'):
        """多重检验校正"""
        p_values = [r['p_value'] for r in results]
        
        if method == 'bonferroni':
            adjusted = [min(p * len(p_values), 1.0) for p in p_values]
        elif method == 'fdr_bh':
            # Benjamini-Hochberg
            adjusted = self.fdr_correction(p_values)
        else:
            adjusted = p_values
        
        for i, r in enumerate(results):
            r['adjusted_pvalue'] = adjusted[i]
            r['log10_pvalue'] = -math.log10(r['p_value']) if r['p_value'] > 0 else 100
            r['log10_adj_pvalue'] = -math.log10(r['adjusted_pvalue']) if r['adjusted_pvalue'] > 0 else 100
        
        return results
    
    def fdr_correction(self, p_values):
        """FDR校正"""
        p_values = np.array(p_values)
        n = len(p_values)
        order = np.argsort(p_values)
        p_sorted = p_values[order]
        
        # 计算校正后的p值
        p_adjusted = p_sorted * n / (np.arange(1, n + 1))
        
        # 确保单调性
        for i in range(n-2, -1, -1):
            if p_adjusted[i] > p_adjusted[i+1]:
                p_adjusted[i] = p_adjusted[i+1]
        
        # 恢复原始顺序
        p_adjusted_original = np.zeros(n)
        p_adjusted_original[order] = p_adjusted
        
        return [min(p, 1.0) for p in p_adjusted_original]
    
    def calculate_confidence_interval(self, k, K, n, N, confidence=0.95):
        """计算置信区间"""
        if k == 0 or n == 0 or K == 0 or N == 0:
            return (0, 0)
        
        # 计算富集比
        enrichment_ratio = (k / n) / (K / N)
        
        # 简化版本的置信区间
        import math
        se = math.sqrt((1/k) + (1/n) + (1/K) + (1/N))
        z_score = 1.96  # 95%置信度
        
        lower = enrichment_ratio * math.exp(-z_score * se)
        upper = enrichment_ratio * math.exp(z_score * se)
        
        return (lower, upper)