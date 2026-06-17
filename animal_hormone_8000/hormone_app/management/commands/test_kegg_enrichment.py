# management/commands/test_kegg_enrichment.py
from django.core.management.base import BaseCommand
from hormone_app.services.kegg_analysis import KeggAnalysisService

class Command(BaseCommand):
    help = '测试KEGG富集分析'
    
    def handle(self, *args, **options):
        # 测试基因列表（细胞周期相关）
        test_genes = [
            "TP53", "CDK1", "CCNB1", "CDKN1A", "RB1", 
            "E2F1", "CCNE1", "CDC25A", "MDM2", "CCND1"
        ]
        
        self.stdout.write(f"测试基因: {test_genes}")
        self.stdout.write("=" * 80)
        
        try:
            analyzer = KeggAnalysisService(organism='hsa')
            results = analyzer.perform_enrichment_analysis(test_genes)
            
            if 'error' in results:
                self.stdout.write(self.style.ERROR(f"错误: {results['error']}"))
                return
            
            # 输出分析结果
            self.stdout.write(self.style.SUCCESS(f"\n✓ 成功映射 {results['mapped_gene_count']} 个基因"))
            self.stdout.write(f"  映射率: {results['mapping_rate']*100:.2f}%")
            self.stdout.write(f"  分析通路: {results['total_pathways_analyzed']}")
            self.stdout.write(f"  显著通路 (p<0.05): {results['significant_pathways']}")
            
            if results['results']:
                self.stdout.write("\n" + "=" * 80)
                self.stdout.write("富集分析结果 (前10个通路):")
                self.stdout.write("=" * 80)
                
                for i, pathway in enumerate(results['results'][:10], 1):
                    self.stdout.write(f"\n{i:2d}. {pathway['pathway_name']}")
                    self.stdout.write(f"     ID: {pathway['pathway_id']}")
                    self.stdout.write(f"     分类: {pathway['category']}")
                    self.stdout.write(f"     P-value: {pathway['p_value']:.2e}")
                    self.stdout.write(f"     校正P值: {pathway['adjusted_pvalue']:.2e}")
                    self.stdout.write(f"     富集分数: {pathway['enrichment_score']:.2f}")
                    self.stdout.write(f"     基因数: {pathway['overlap_count']}/{pathway['pathway_gene_count']}")
                    self.stdout.write(f"     基因列表: {', '.join(pathway['overlap_genes'][:5])}")
                
                # 期望看到细胞周期通路显著富集
                cell_cycle = next((p for p in results['results'] if 'Cell cycle' in p['pathway_name']), None)
                if cell_cycle:
                    self.stdout.write(self.style.SUCCESS(
                        f"\n✓ 细胞周期通路显著富集! (p={cell_cycle['adjusted_pvalue']:.2e})"
                    ))
                
        except Exception as e:
            self.stdout.write(self.style.ERROR(f"错误: {e}"))
            import traceback
            traceback.print_exc()