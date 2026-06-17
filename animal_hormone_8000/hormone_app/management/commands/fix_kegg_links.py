# management/commands/fix_kegg_links.py
import json
from pathlib import Path
from django.core.management.base import BaseCommand
from django.db import transaction
from hormone_app.models import KeggGene, KeggPathway, KeggGenePathway, KeggOrganism

class Command(BaseCommand):
    help = '修复KEGG基因-通路关联'
    
    def add_arguments(self, parser):
        parser.add_argument('--organism', type=str, default='hsa')
        parser.add_argument('--data-dir', type=str, default='data/kegg')
    
    def handle(self, *args, **options):
        organism = options['organism']
        data_dir = Path(options['data_dir'])
        
        try:
            org = KeggOrganism.objects.get(code=organism)
            
            # 1. 首先检查有哪些关联文件
            self.stdout.write("检查关联文件...")
            
            possible_files = [
                data_dir / f"relations_{organism}.json",
                data_dir / f"links_{organism}.json",
                data_dir / f"relations_kgml_{organism}.json",
                data_dir / "relations_hsa.json",
                data_dir / "links_hsa.json"
            ]
            
            relation_file = None
            for file_path in possible_files:
                if file_path.exists():
                    relation_file = file_path
                    self.stdout.write(f"找到关联文件: {file_path}")
                    break
            
            if not relation_file:
                self.stdout.write(self.style.ERROR("未找到关联文件，需要重新下载"))
                self.download_links(organism, data_dir)
                return
            
            # 2. 读取并分析关联文件
            with open(relation_file, 'r') as f:
                relations = json.load(f)
            
            self.stdout.write(f"\n关联文件分析:")
            self.stdout.write(f"  总记录数: {len(relations)}")
            
            if len(relations) == 0:
                self.stdout.write(self.style.ERROR("关联文件为空"))
                self.download_links(organism, data_dir)
                return
            
            # 显示前几条记录的结构
            self.stdout.write("\n  前5条记录结构:")
            for i, rel in enumerate(relations[:5]):
                self.stdout.write(f"    {i+1}. {rel}")
            
            # 3. 获取所有基因和通路
            genes = {g.gene_id: g for g in KeggGene.objects.filter(organism=org)}
            pathways = {p.pathway_id: p for p in KeggPathway.objects.all()}
            
            self.stdout.write(f"\n数据库状态:")
            self.stdout.write(f"  基因数: {len(genes)}")
            self.stdout.write(f"  通路数: {len(pathways)}")
            
            # 显示一些基因示例
            self.stdout.write("\n  基因示例 (前5个):")
            for gene_id in list(genes.keys())[:5]:
                self.stdout.write(f"    {gene_id}: {genes[gene_id].symbol}")
            
            # 显示一些通路示例
            self.stdout.write("\n  通路示例 (前5个):")
            for pathway_id in list(pathways.keys())[:5]:
                self.stdout.write(f"    {pathway_id}: {pathways[pathway_id].name[:30]}...")
            
            # 4. 尝试匹配关联
            matched = 0
            format_issues = []
            
            for rel in relations[:1000]:  # 只测试前1000条
                # 尝试不同的字段名
                gene_id = None
                pathway_id = None
                
                if isinstance(rel, dict):
                    # 尝试各种可能的字段名
                    if 'gene_id' in rel:
                        gene_id = str(rel['gene_id'])
                    elif 'gene' in rel:
                        gene_id = str(rel['gene'])
                    elif 'gene_symbol' in rel:
                        # 可能需要通过符号查找基因ID
                        symbol = rel['gene_symbol']
                        gene = KeggGene.objects.filter(organism=org, symbol=symbol).first()
                        if gene:
                            gene_id = gene.gene_id
                    
                    if 'pathway_id' in rel:
                        pathway_id = rel['pathway_id']
                    elif 'pathway' in rel:
                        pathway_id = rel['pathway']
                
                if gene_id and pathway_id:
                    if gene_id in genes and pathway_id in pathways:
                        matched += 1
                    else:
                        format_issues.append({
                            'gene_id': gene_id,
                            'pathway_id': pathway_id,
                            'gene_exists': gene_id in genes,
                            'pathway_exists': pathway_id in pathways
                        })
            
            self.stdout.write(f"\n匹配测试 (前1000条):")
            self.stdout.write(f"  可匹配: {matched}")
            self.stdout.write(f"  格式问题: {len(format_issues)}")
            
            if format_issues:
                self.stdout.write("\n  格式问题示例:")
                for issue in format_issues[:5]:
                    self.stdout.write(f"    Gene: {issue['gene_id']} (存在: {issue['gene_exists']}), "
                                    f"Pathway: {issue['pathway_id']} (存在: {issue['pathway_exists']})")
            
            # 5. 重新下载关联数据
            if matched == 0:
                self.stdout.write(self.style.WARNING("\n无法匹配任何关联，重新下载..."))
                self.download_links(organism, data_dir)
                
        except KeggOrganism.DoesNotExist:
            self.stdout.write(self.style.ERROR(f"生物体 {organism} 不存在"))
    
    def download_links(self, organism, data_dir):
        """重新下载关联数据"""
        import requests
        
        self.stdout.write("\n重新下载基因-通路关联...")
        
        url = f"http://rest.kegg.jp/link/pathway/{organism}"
        
        try:
            response = requests.get(url, timeout=60)
            
            if response.status_code != 200:
                self.stdout.write(self.style.ERROR(f"下载失败: {response.status_code}"))
                return
            
            content = response.text.strip()
            lines = content.split('\n')
            
            self.stdout.write(f"下载了 {len(lines)} 行数据")
            
            relations = []
            for line in lines:
                if line:
                    parts = line.split('\t')
                    if len(parts) == 2:
                        gene_part, pathway_part = parts
                        
                        # 提取基因ID (格式: hsa:1)
                        if ':' in gene_part:
                            gene_id = gene_part.split(':')[1]
                        else:
                            gene_id = gene_part
                        
                        # 通路ID保持原样
                        pathway_id = pathway_part
                        
                        relations.append({
                            'gene_id': gene_id,
                            'pathway_id': pathway_id
                        })
            
            # 保存为新文件
            new_file = data_dir / f"links_fixed_{organism}.json"
            with open(new_file, 'w') as f:
                json.dump(relations, f, indent=2)
            
            self.stdout.write(self.style.SUCCESS(f"已保存修复后的关联文件: {new_file}"))
            self.stdout.write(f"共 {len(relations)} 条关联")
            
            # 立即导入
            self.import_links(organism, relations)
            
        except Exception as e:
            self.stdout.write(self.style.ERROR(f"下载失败: {e}"))
    
    def import_links(self, organism, relations):
        """导入关联"""
        try:
            org = KeggOrganism.objects.get(code=organism)
            genes = {g.gene_id: g for g in KeggGene.objects.filter(organism=org)}
            pathways = {p.pathway_id: p for p in KeggPathway.objects.all()}
            
            self.stdout.write(f"\n开始导入关联...")
            
            with transaction.atomic():
                # 删除旧关联
                deleted = KeggGenePathway.objects.filter(gene__organism=org).delete()[0]
                self.stdout.write(f"删除 {deleted} 条旧关联")
                
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
                
        except Exception as e:
            self.stdout.write(self.style.ERROR(f"导入失败: {e}"))
    
    def update_counts(self, org):
        """更新计数"""
        from django.db import connection
        
        self.stdout.write("更新计数...")
        
        with connection.cursor() as cursor:
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