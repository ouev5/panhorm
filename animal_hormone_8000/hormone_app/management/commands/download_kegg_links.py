# management/commands/download_kegg_links.py
import os
import requests
import json
from pathlib import Path
from django.core.management.base import BaseCommand
from django.db import transaction, connection
from tqdm import tqdm
from ...models import KeggOrganism, KeggPathway, KeggGene, KeggGenePathway

class Command(BaseCommand):
    help = '下载KEGG基因-通路关联链接'
    
    def add_arguments(self, parser):
        parser.add_argument('--organism', type=str, default='hsa', 
                          help='要下载的生物体代码')
    
    def handle(self, *args, **options):
        organism = options['organism']
        
        self.stdout.write(f"下载KEGG基因-通路关联，生物体: {organism}")
        
        try:
            org_obj = KeggOrganism.objects.get(code=organism)
            
            # 1. 下载关联数据
            relations = self.download_link_data(organism)
            
            # 2. 导入到数据库
            self.import_links_to_db(organism, org_obj, relations)
            
        except KeggOrganism.DoesNotExist:
            self.stdout.write(self.style.ERROR(f"生物体 {organism} 不存在"))
            return
    
    def download_link_data(self, organism):
        """下载链接数据"""
        self.stdout.write("下载基因-通路关联链接...")
        
        url = f"http://rest.kegg.jp/link/{organism}/pathway"
        
        try:
            response = requests.get(url, timeout=60)
            
            if response.status_code != 200:
                self.stdout.write(f"错误: API返回状态码 {response.status_code}")
                return []
            
            content = response.text.strip()
            if not content:
                self.stdout.write("警告: API返回空数据")
                return []
            
            # 解析数据
            relations = []
            lines = content.split('\n')
            
            self.stdout.write(f"解析 {len(lines)} 行数据...")
            
            for line in tqdm(lines, desc="解析关联"):
                if line:
                    parts = line.split('\t')
                    if len(parts) == 2:
                        gene_part, pathway_part = parts
                        
                        # 提取基因ID (格式: hsa:1)
                        if ':' in gene_part:
                            gene_id = gene_part.split(':')[1]  # 取冒号后的部分
                        else:
                            gene_id = gene_part.replace(f'{organism}:', '')
                        
                        # 提取通路ID (格式: path:hsa00010)
                        if ':' in pathway_part:
                            pathway_id = pathway_part  # 保持原样
                        else:
                            pathway_id = f'path:{pathway_part}' if not pathway_part.startswith('path:') else pathway_part
                        
                        relations.append({
                            'gene_id': gene_id,
                            'pathway_id': pathway_id
                        })
            
            # 保存到文件
            data_dir = Path('data/kegg')
            data_dir.mkdir(parents=True, exist_ok=True)
            
            relation_file = data_dir / f"links_{organism}.json"
            with open(relation_file, 'w') as f:
                json.dump(relations, f, indent=2)
            
            self.stdout.write(f"找到 {len(relations)} 条关联")
            return relations
            
        except Exception as e:
            self.stdout.write(f"错误: 下载关联数据失败: {e}")
            return []
    
    def import_links_to_db(self, organism, org_obj, relations):
        """导入链接到数据库"""
        if not relations:
            self.stdout.write("没有关联数据需要导入")
            return
        
        self.stdout.write(f"导入 {len(relations)} 条关联到数据库...")
        
        try:
            with transaction.atomic():
                # 获取所有基因ID和通路ID
                gene_ids = list(set([r['gene_id'] for r in relations]))
                pathway_ids = list(set([r['pathway_id'] for r in relations]))
                
                self.stdout.write(f"去重后: {len(gene_ids)} 个基因, {len(pathway_ids)} 条通路")
                
                # 批量获取基因和通路对象
                genes_dict = {gene.gene_id: gene for gene in 
                            KeggGene.objects.filter(gene_id__in=gene_ids, organism=org_obj)}
                pathways_dict = {pathway.pathway_id: pathway for pathway in 
                               KeggPathway.objects.filter(pathway_id__in=pathway_ids)}
                
                self.stdout.write(f"数据库匹配: {len(genes_dict)} 个基因, {len(pathways_dict)} 条通路")
                
                # 准备批量插入
                relations_to_create = []
                skipped_missing = 0
                
                for rel in tqdm(relations, desc="准备关联数据"):
                    gene_id = rel['gene_id']
                    pathway_id = rel['pathway_id']
                    
                    if gene_id in genes_dict and pathway_id in pathways_dict:
                        relations_to_create.append(
                            KeggGenePathway(
                                gene=genes_dict[gene_id],
                                pathway=pathways_dict[pathway_id],
                                relation_type='member',
                                evidence='KEGG'
                            )
                        )
                    else:
                        skipped_missing += 1
                
                if skipped_missing > 0:
                    self.stdout.write(f"跳过 {skipped_missing} 条关联（基因或通路不存在）")
                
                # 批量创建关联（使用ignore_conflicts避免重复）
                if relations_to_create:
                    self.stdout.write(f"批量创建 {len(relations_to_create)} 条关联...")
                    
                    # 先删除可能存在的旧关联
                    KeggGenePathway.objects.filter(gene__organism=org_obj).delete()
                    
                    # 批量创建新关联
                    batch_size = 1000
                    for i in range(0, len(relations_to_create), batch_size):
                        batch = relations_to_create[i:i+batch_size]
                        KeggGenePathway.objects.bulk_create(batch, ignore_conflicts=True)
                
                # 更新计数
                self.update_counts(org_obj)
                
                # 显示最终统计
                total_relations = KeggGenePathway.objects.filter(gene__organism=org_obj).count()
                self.stdout.write(self.style.SUCCESS(
                    f"\n导入完成:\n"
                    f"  总关联数: {total_relations}\n"
                    f"  基因平均通路数: {total_relations / len(genes_dict) if genes_dict else 0:.1f}\n"
                    f"  通路平均基因数: {total_relations / len(pathways_dict) if pathways_dict else 0:.1f}\n"
                ))
                
        except Exception as e:
            self.stdout.write(self.style.ERROR(f"导入失败: {e}"))
            import traceback
            traceback.print_exc()
    
    def update_counts(self, org_obj):
        """更新计数"""
        self.stdout.write("更新通路和基因的关联计数...")
        
        # 使用原生SQL更新通路基因计数
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
            """, [org_obj.code])
            
            # 更新基因通路计数
            cursor.execute("""
                UPDATE hormone_app_kegggene g
                SET pathway_count = (
                    SELECT COUNT(DISTINCT gp.pathway_id)
                    FROM hormone_app_kegggenepathway gp
                    WHERE gp.gene_id = g.gene_id
                )
                WHERE g.organism_id = %s
            """, [org_obj.code])