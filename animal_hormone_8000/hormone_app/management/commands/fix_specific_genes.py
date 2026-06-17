# management/commands/fix_specific_genes.py
from django.core.management.base import BaseCommand
from django.db import transaction
from hormone_app.models import KeggGene, KeggOrganism

class Command(BaseCommand):
    help = '修复特定基因的符号'
    
    def handle(self, *args, **options):
        try:
            org = KeggOrganism.objects.get(code='hsa')
            
            # 需要修复的基因映射表
            fix_map = {
                # 细胞周期相关基因
                '7157': 'TP53',      # TP53
                '983': 'CDK1',       # CDK1
                '891': 'CCNB1',      # CCNB1
                '1026': 'CDKN1A',    # CDKN1A
                '5925': 'RB1',       # RB1 (之前是5542，但5925也是RB1)
                '1869': 'E2F1',      # E2F1
                '898': 'CCNE1',      # CCNE1
                '993': 'CDC25A',     # CDC25A
                '4193': 'MDM2',      # MDM2
                '595': 'CCND1',      # CCND1
                
                # 其他常见基因
                '672': 'BRCA1',      # BRCA1
                '675': 'BRCA2',      # BRCA2
                '1956': 'EGFR',      # EGFR
                '3845': 'KRAS',      # KRAS
                '5290': 'PIK3CA',    # PIK3CA
                '207': 'AKT1',       # AKT1
                '5728': 'PTEN',      # PTEN
                '4609': 'MYC',       # MYC
                '1019': 'CDK4',      # CDK4
                '1021': 'CDK6',      # CDK6
            }
            
            self.stdout.write("修复特定基因符号...")
            
            with transaction.atomic():
                updated = 0
                not_found = []
                already_correct = 0
                
                for gene_id, correct_symbol in fix_map.items():
                    try:
                        gene = KeggGene.objects.get(gene_id=gene_id, organism=org)
                        old_symbol = gene.symbol
                        
                        # 如果当前符号已经是GENE_ID格式或者不是正确符号
                        if old_symbol != correct_symbol:
                            gene.symbol = correct_symbol
                            gene.save()
                            updated += 1
                            self.stdout.write(f"  ✓ {gene_id}: {old_symbol} -> {correct_symbol}")
                        else:
                            already_correct += 1
                            self.stdout.write(f"  - {gene_id}: 已经是 {correct_symbol}")
                            
                    except KeggGene.DoesNotExist:
                        not_found.append(gene_id)
                
                self.stdout.write(self.style.SUCCESS(f"\n成功修复 {updated} 个基因"))
                self.stdout.write(f"已正确 {already_correct} 个基因")
                if not_found:
                    self.stdout.write(self.style.WARNING(f"未找到的基因ID: {not_found}"))
                
        except KeggOrganism.DoesNotExist:
            self.stdout.write(self.style.ERROR("生物体hsa不存在"))