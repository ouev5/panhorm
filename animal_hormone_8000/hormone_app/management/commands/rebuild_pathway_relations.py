# management/commands/rebuild_pathway_relations.py
from django.core.management.base import BaseCommand
from django.db import transaction, connection
from hormone_app.models import KeggOrganism, KeggPathway, KeggGene, KeggGenePathway
import json
from pathlib import Path

class Command(BaseCommand):
    help = '重建通路基因关联'
    
    def add_arguments(self, parser):
        parser.add_argument('--organism', type=str, default='hsa')
        parser.add_argument('--data-dir', type=str, default='data/kegg')
    
    def handle(self, *args, **options):
        organism = options['organism']
        data_dir = Path(options['data_dir'])
        
        try:
            org = KeggOrganism.objects.get(code=organism)
            
            self.stdout.write(f"重建通路基因关联...")
            
            # 读取关联数据
            relation_file = data_dir / f"links_{organism}.json"
            if not relation_file.exists():
                self.stdout.write(self.style.ERROR(f"关联文件不存在: {relation_file}"))
                return
            
            with open(relation_file, 'r') as f:
                relations = json.load(f)
            
            self.stdout.write(f"加载了 {len(relations)} 条关联")
            
            with transaction.atomic():
                # 删除旧的关联
                deleted = KeggGenePathway.objects.filter(gene__organism=org).delete()[0]
                self.stdout.write(f"删除 {deleted} 条旧关联")
                
                # 获取所有基因和通路
                genes = {g.gene_id: g for g in KeggGene.objects.filter(organism=org)}
                pathways = {p.pathway_id: p for p in KeggPathway.objects.all()}
                
                self.stdout.write(f"数据库中有 {len(genes)} 个基因, {len(pathways)} 条通路")
                
                # 创建新关联
                created = 0
                skipped = 0
                
                for rel in relations:
                    gene_id = rel['gene_id']
                    pathway_id = rel['pathway_id']
                    
                    if gene_id in genes and pathway_id in pathways:
                        KeggGenePathway.objects.create(
                            gene=genes[gene_id],
                            pathway=pathways[pathway_id],
                            relation_type='member',
                            evidence='KEGG'
                        )
                        created += 1
                    else:
                        skipped += 1
                    
                    if created % 1000 == 0:
                        self.stdout.write(f"  已创建 {created} 条关联...")
                
                self.stdout.write(self.style.SUCCESS(f"\n成功创建 {created} 条关联"))
                self.stdout.write(f"跳过 {skipped} 条关联")
                
                # 更新计数
                self.update_counts(org)
                
        except KeggOrganism.DoesNotExist:
            self.stdout.write(self.style.ERROR(f"生物体 {organism} 不存在"))
    
    def update_counts(self, org):
        """更新通路和基因的计数"""
        self.stdout.write("更新计数...")
        
        with connection.cursor() as cursor:
            # 更新通路基因计数
            cursor.execute("""
                UPDATE hormone_app_keggpathway p
                SET gene_count = (
                    SELECT COUNT(DISTINCT g.gene_id)
                    FROM hormone_app_kegggenepathway gp
                    JOIN hormone_app_kegggene g ON gp.gene_id = g.gene_id
                    WHERE gp.pathway_id = p.pathway_id
                    AND g.organism_id = %s
                )
            """, [org.code])
            
            # 更新基因通路计数
            cursor.execute("""
                UPDATE hormone_app_kegggene g
                SET pathway_count = (
                    SELECT COUNT(DISTINCT gp.pathway_id)
                    FROM hormone_app_kegggenepathway gp
                    WHERE gp.gene_id = g.gene_id
                )
                WHERE g.organism_id = %s
            """, [org.code])
        
        self.stdout.write("计数更新完成")