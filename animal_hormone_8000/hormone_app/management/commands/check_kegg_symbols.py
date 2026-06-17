# management/commands/check_kegg_symbols.py
from django.core.management.base import BaseCommand
from hormone_app.models import KeggGene, KeggOrganism

class Command(BaseCommand):
    help = '检查KEGG基因符号格式'
    
    def handle(self, *args, **options):
        try:
            org = KeggOrganism.objects.get(code='hsa')
            
            # 检查一些基因样本
            sample_genes = KeggGene.objects.filter(organism=org)[:20]
            
            self.stdout.write("=" * 80)
            self.stdout.write("数据库中基因符号示例:")
            self.stdout.write("=" * 80)
            
            for gene in sample_genes:
                self.stdout.write(f"ID: {gene.gene_id:10} | Symbol: '{gene.symbol}'")
            
            # 检查特定基因
            test_genes = ["TP53", "CDK1", "CCNB1", "CDKN1A", "RB1", "E2F1", "CCNE1", "CDC25A", "MDM2", "CCND1"]
            
            self.stdout.write("\n" + "=" * 80)
            self.stdout.write("测试基因匹配情况:")
            self.stdout.write("=" * 80)
            
            for gene_symbol in test_genes:
                # 精确匹配
                exact = KeggGene.objects.filter(organism=org, symbol=gene_symbol)
                # 模糊匹配
                fuzzy = KeggGene.objects.filter(organism=org, symbol__icontains=gene_symbol)
                # 包含匹配
                contains = KeggGene.objects.filter(organism=org, symbol__contains=gene_symbol)
                
                self.stdout.write(f"\n{gene_symbol}:")
                self.stdout.write(f"  - 精确匹配: {exact.count()}")
                if exact.count() > 0:
                    for g in exact:
                        self.stdout.write(f"    * {g.gene_id}: '{g.symbol}'")
                
                self.stdout.write(f"  - 模糊匹配: {fuzzy.count()}")
                if fuzzy.count() > 0 and exact.count() == 0:
                    for g in fuzzy[:3]:
                        self.stdout.write(f"    * {g.gene_id}: '{g.symbol}'")
                
                self.stdout.write(f"  - 包含匹配: {contains.count()}")
            
            # 统计
            total_genes = KeggGene.objects.filter(organism=org).count()
            null_symbols = KeggGene.objects.filter(organism=org, symbol__isnull=True).count()
            empty_symbols = KeggGene.objects.filter(organism=org, symbol='').count()
            
            self.stdout.write("\n" + "=" * 80)
            self.stdout.write("统计信息:")
            self.stdout.write("=" * 80)
            self.stdout.write(f"总基因数: {total_genes}")
            self.stdout.write(f"空符号: {null_symbols + empty_symbols}")
            
        except KeggOrganism.DoesNotExist:
            self.stdout.write(self.style.ERROR("生物体hsa不存在"))