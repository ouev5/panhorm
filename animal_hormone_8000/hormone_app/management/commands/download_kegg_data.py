# management/commands/download_kegg_data.py
import os
import requests
import gzip
import xml.etree.ElementTree as ET
from django.core.management.base import BaseCommand
from django.db import transaction
from django.core.exceptions import ValidationError
from tqdm import tqdm
import io
import json
from pathlib import Path
import time
import re
# 注意：根据实际项目结构调整导入路径
from ...models import KeggOrganism, KeggPathway, KeggGene, KeggGenePathway

class Command(BaseCommand):
    help = '下载并导入KEGG数据库到MySQL'
    
    def add_arguments(self, parser):
        parser.add_argument('--organism', type=str, default='hsa', 
                          help='要下载的生物体代码，如hsa, mmu, rno等')
        parser.add_argument('--update', action='store_true', 
                          help='更新现有数据而不是重新创建')
        parser.add_argument('--data-dir', type=str, default='data/kegg',
                          help='数据存储目录')
        parser.add_argument('--skip-relations', action='store_true',
                          help='跳过基因-通路关联下载')
    
    def handle(self, *args, **options):
        self.data_dir = Path(options['data_dir'])
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
        organism_code = options['organism']
        self.update_mode = options['update']
        self.skip_relations = options['skip_relations']
        
        self.stdout.write(f"开始下载KEGG数据，生物体: {organism_code}")
        
        # 0. 先创建或获取生物体对象
        org_obj = self.get_or_create_organism(organism_code)
        
        # 1. 下载通路列表
        pathways = self.download_pathway_list(organism_code)
        
        # 2. 下载基因列表
        genes = self.download_gene_list(organism_code)
        
        # 3. 下载基因-通路关联
        if not self.skip_relations:
            relations = self.download_gene_pathway_relations(organism_code, pathways, genes)
        else:
            self.stdout.write("跳过基因-通路关联下载")
            relations = []
        
        # 4. 导入数据库
        self.import_to_database(organism_code, org_obj, pathways, genes)
        
        self.stdout.write(self.style.SUCCESS('KEGG数据导入完成！'))
    
    def get_or_create_organism(self, organism_code):
        """创建或获取生物体对象"""
        self.stdout.write(f"创建/获取生物体对象: {organism_code}")
        
        try:
            # 尝试获取现有对象
            org_obj = KeggOrganism.objects.get(code=organism_code)
            self.stdout.write(f"使用现有生物体对象: {org_obj.name}")
        except KeggOrganism.DoesNotExist:
            # 创建新的生物体对象 - 使用正确的字段
            org_obj = KeggOrganism.objects.create(
                code=organism_code,
                name=self.get_organism_name(organism_code),
                scientific_name=self.get_scientific_name(organism_code),
                description=f"{self.get_organism_name(organism_code)} ({self.get_scientific_name(organism_code)})"
            )
            self.stdout.write(f"创建新生物体对象: {org_obj.name}")
        
        return org_obj
    
    def download_pathway_list(self, organism):
        """下载通路列表"""
        self.stdout.write("下载通路列表...")
        
        url = f"http://rest.kegg.jp/list/pathway/{organism}"
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            pathways = []
            for line in response.text.strip().split('\n'):
                if line:
                    parts = line.split('\t', 1)
                    if len(parts) == 2:
                        pathway_id, name = parts
                        # 为每个通路添加分类
                        category = self.categorize_pathway(pathway_id, name)
                        pathways.append({
                            'id': pathway_id,
                            'name': name,
                            'organism': organism,
                            'category': category
                        })
            
            # 保存到文件
            pathway_file = self.data_dir / f"pathways_{organism}.json"
            with open(pathway_file, 'w') as f:
                json.dump(pathways, f, indent=2)
            
            self.stdout.write(f"找到 {len(pathways)} 条通路")
            return pathways
        except Exception as e:
            self.stdout.write(f"错误: 下载通路列表失败: {e}")
            return []
    
    def download_gene_list(self, organism):
        """下载基因列表"""
        self.stdout.write("下载基因列表...")
        
        # 方法1: 从KEGG API下载
        url = f"http://rest.kegg.jp/list/{organism}"
        try:
            response = requests.get(url, timeout=60)
            response.raise_for_status()
            
            genes = []
            total_lines = response.text.strip().split('\n')
            
            for line in tqdm(total_lines, desc="解析基因"):
                if line:
                    parts = line.split('\t', 1)
                    if len(parts) >= 2:
                        gene_id = parts[0]  # 格式如: hsa:1
                        description = parts[1]
                        
                        # 提取基因符号和名称
                        name_parts = description.split(';')
                        symbol = name_parts[0].strip() if name_parts else gene_id
                        
                        # 清理符号 - 移除制表符等特殊字符
                        symbol = self.clean_gene_symbol(symbol)
                        
                        # 提取Entrez ID等信息
                        entrez_id = ''
                        uniprot_id = ''
                        ensembl_id = ''
                        for part in name_parts:
                            part = part.strip()
                            if 'Entrez:' in part:
                                entrez_id = part.split(':')[-1].strip()
                            elif 'UniProt:' in part:
                                uniprot_id = part.split(':')[-1].strip()
                            elif 'Ensembl:' in part:
                                ensembl_id = part.split(':')[-1].strip()
                        
                        # 只保留基因ID的最后部分（如hsa:1 -> 1）
                        gene_id_short = gene_id.split(':')[-1] if ':' in gene_id else gene_id
                        
                        # 如果符号仍然有问题，使用基因ID作为符号
                        if not symbol or symbol == '' or '\t' in symbol:
                            symbol = f"GENE_{gene_id_short}"
                        
                        # 截断过长的symbol
                        if len(symbol) > 200:
                            symbol = symbol[:197] + "..."
                        
                        genes.append({
                            'gene_id': gene_id_short,  # 使用短ID
                            'kegg_id': gene_id,  # 保存完整KEGG ID
                            'symbol': symbol,
                            'name': description,
                            'entrez_id': entrez_id,
                            'uniprot_id': uniprot_id,
                            'ensembl_id': ensembl_id,
                            'organism': organism
                        })
            
            # 保存到文件
            gene_file = self.data_dir / f"genes_{organism}.json"
            with open(gene_file, 'w') as f:
                json.dump(genes, f, indent=2)
            
            self.stdout.write(f"找到 {len(genes)} 个基因")
            return genes
        except Exception as e:
            self.stdout.write(f"错误: 下载基因列表失败: {e}")
            return []
    
    def clean_gene_symbol(self, symbol):
        """清理基因符号"""
        if not symbol:
            return symbol
            
        # 移除制表符、换行符等
        symbol = symbol.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ')
        
        # 移除多余空格
        symbol = re.sub(r'\s+', ' ', symbol).strip()
        
        # 如果包含类似"CDS\tY:complement"的结构，提取有意义的部分
        if ':' in symbol:
            # 提取冒号前的内容
            parts = symbol.split(':')
            if parts[0]:
                symbol = parts[0].strip()
        
        # 如果包含"CDS"、"tRNA"等前缀，尝试提取后面的内容
        patterns = [
            r'CDS\s+(.+)',
            r'tRNA\s+(.+)',
            r'ncRNA\s+(.+)',
            r'rRNA\s+(.+)'
        ]
        
        for pattern in patterns:
            match = re.match(pattern, symbol)
            if match:
                extracted = match.group(1).strip()
                if extracted and extracted != 'Y' and extracted != 'MT':
                    symbol = extracted
                    break
        
        # 如果符号仍然包含不必要的内容，尝试简化
        if symbol in ['Y', 'MT', 'complement'] or len(symbol) < 2:
            return f"GENE_{hash(symbol) % 10000}"  # 生成一个唯一的符号
        
        return symbol[:200]  # 确保不超过200字符
    
    def download_gene_pathway_relations(self, organism, pathways, genes):
        """下载基因-通路关联 - 使用更可靠的方法"""
        self.stdout.write("下载基因-通路关联...")
        
        relations = []
        gene_dict = {g['gene_id']: g for g in genes}
        
        try:
            # 方法1：逐个通路下载关联
            for pathway in tqdm(pathways, desc="获取通路基因"):
                pathway_id_short = pathway['id'].replace('path:', '')
                
                # 使用更稳定的API
                url = f"http://rest.kegg.jp/get/{pathway_id_short}"
                try:
                    response = requests.get(url, timeout=30)
                    if response.status_code == 200:
                        # 解析返回的文本，提取基因信息
                        lines = response.text.split('\n')
                        gene_section = False
                        
                        for line in lines:
                            if line.startswith('GENE'):
                                gene_section = True
                            elif gene_section and line.startswith(' '):
                                # 解析基因行，格式如:       1       hsa:1        A1BG
                                parts = line.strip().split()
                                if len(parts) >= 3:
                                    # 提取基因ID（如hsa:1）
                                    possible_gene = parts[1]
                                    if ':' in possible_gene and possible_gene.startswith(f'{organism}:'):
                                        gene_id = possible_gene.split(':')[-1]
                                        
                                        if gene_id in gene_dict:
                                            relations.append({
                                                'gene_id': gene_id,
                                                'gene_symbol': gene_dict[gene_id]['symbol'],
                                                'pathway_id': pathway['id'],
                                                'pathway_name': pathway['name']
                                            })
                            elif gene_section and not line.startswith(' '):
                                # 结束基因部分
                                break
                except Exception as e:
                    self.stdout.write(f"警告: 获取通路 {pathway_id_short} 失败: {e}")
                    continue
            
            if len(relations) == 0:
                # 方法2：使用link API作为备用
                self.stdout.write("尝试备用方法获取关联...")
                url = f"http://rest.kegg.jp/link/{organism}/pathway"
                try:
                    response = requests.get(url, timeout=60)
                    if response.status_code == 200:
                        for line in response.text.strip().split('\n'):
                            if line:
                                parts = line.split('\t')
                                if len(parts) == 2:
                                    gene_part, pathway_part = parts
                                    
                                    # 提取基因ID
                                    if ':' in gene_part:
                                        gene_id = gene_part.split(':')[-1]
                                    else:
                                        gene_id = gene_part.replace(f'{organism}:', '')
                                    
                                    # 提取通路ID
                                    if ':' in pathway_part:
                                        pathway_id = pathway_part.split(':')[-1]
                                    else:
                                        pathway_id = pathway_part.replace('path:', '')
                                    
                                    full_pathway_id = f'path:{pathway_id}'
                                    
                                    if gene_id in gene_dict and full_pathway_id in [p['id'] for p in pathways]:
                                        relations.append({
                                            'gene_id': gene_id,
                                            'gene_symbol': gene_dict[gene_id]['symbol'],
                                            'pathway_id': full_pathway_id,
                                            'pathway_name': next((p['name'] for p in pathways if p['id'] == full_pathway_id), '')
                                        })
                except Exception as e:
                    self.stdout.write(f"警告: 备用方法也失败: {e}")
            
            # 保存到文件
            relation_file = self.data_dir / f"relations_{organism}.json"
            with open(relation_file, 'w') as f:
                json.dump(relations, f, indent=2)
            
            self.stdout.write(f"找到 {len(relations)} 条关联")
            return relations
            
        except Exception as e:
            self.stdout.write(f"错误: 下载基因-通路关联失败: {e}")
            return []
    
    def import_to_database(self, organism, org_obj, pathways, genes):
        """导入数据到MySQL数据库 - 修复基因重复问题"""
        self.stdout.write("导入数据到MySQL数据库...")
        
        try:
            with transaction.atomic():
                # 如果是更新模式，先删除旧数据
                if self.update_mode:
                    self.stdout.write("更新模式：删除旧数据...")
                    KeggGenePathway.objects.filter(gene__organism=org_obj).delete()
                    KeggGene.objects.filter(organism=org_obj).delete()
                
                # 1. 导入通路
                pathway_objs = {}
                self.stdout.write(f"导入通路 ({len(pathways)}条)...")
                
                for pathway in tqdm(pathways, desc="导入通路"):
                    pathway_id = pathway['id']
                    try:
                        # 准备通路数据
                        pathway_data = {
                            'name': pathway['name'][:500],
                            'description': '',
                            'category': pathway.get('category', ''),
                            'subcategory': '',
                            'gene_count': 0,
                            'compound_count': 0,
                            'reaction_count': 0,
                            'kegg_url': f"https://www.genome.jp/pathway/{pathway_id}",
                            'image_url': f"https://www.genome.jp/pathway/{pathway_id}/{pathway_id}.png"
                        }
                        
                        obj, created = KeggPathway.objects.get_or_create(
                            pathway_id=pathway_id,
                            defaults=pathway_data
                        )
                        
                        # 如果已存在，更新数据
                        if not created:
                            for key, value in pathway_data.items():
                                if hasattr(obj, key):
                                    setattr(obj, key, value)
                            obj.save()
                        
                        pathway_objs[pathway_id] = obj
                    except Exception as e:
                        self.stdout.write(f"警告: 导入通路 {pathway_id} 失败: {e}")
                        continue
                
                # 2. 导入基因 - 处理重复符号问题
                gene_objs = {}
                symbol_counter = {}  # 记录每个符号出现的次数
                self.stdout.write(f"导入基因 ({len(genes)}条)...")
                
                # 第一遍：统计符号使用情况
                for gene in genes:
                    symbol = gene.get('symbol', '')
                    if symbol:
                        symbol_counter[symbol] = symbol_counter.get(symbol, 0) + 1
                
                # 第二遍：导入基因，处理重复符号
                imported_count = 0
                duplicate_count = 0
                
                for gene in tqdm(genes, desc="导入基因"):
                    try:
                        gene_id = gene['gene_id']
                        symbol = gene.get('symbol', '')
                        
                        # 如果符号重复，添加后缀
                        if symbol and symbol_counter.get(symbol, 0) > 1:
                            # 查找已导入的相同符号的基因
                            existing_count = KeggGene.objects.filter(
                                organism=org_obj, 
                                symbol=symbol
                            ).count()
                            
                            if existing_count > 0:
                                # 添加后缀确保唯一性
                                new_symbol = f"{symbol}_{existing_count + 1}"
                                gene['symbol'] = new_symbol
                                duplicate_count += 1
                        
                        # 准备入库数据
                        gene_data = {
                            'symbol': gene.get('symbol', ''),
                            'name': gene.get('name', ''),
                            'entrez_id': gene.get('entrez_id', ''),
                            'uniprot_id': gene.get('uniprot_id', ''),
                            'ensembl_id': gene.get('ensembl_id', ''),
                            'organism': org_obj,
                            'pathway_count': 0
                        }
                        
                        # 如果symbol为空，使用gene_id
                        if not gene_data['symbol']:
                            gene_data['symbol'] = f"GENE_{gene_id}"
                        
                        # 确保symbol不超过200字符
                        if len(gene_data['symbol']) > 200:
                            gene_data['symbol'] = gene_data['symbol'][:200]
                        
                        # 尝试创建基因记录
                        try:
                            obj, created = KeggGene.objects.get_or_create(
                                gene_id=gene_id,
                                organism=org_obj,
                                defaults=gene_data
                            )
                            
                            # 如果已存在，更新数据（但不更新symbol以避免验证错误）
                            if not created:
                                # 只更新除symbol外的其他字段
                                for key, value in gene_data.items():
                                    if key not in ['symbol', 'organism'] and hasattr(obj, key):
                                        setattr(obj, key, value)
                                obj.save()
                            
                            gene_objs[gene_id] = obj
                            imported_count += 1
                            
                        except ValidationError as ve:
                            # 如果验证失败，尝试使用不同的symbol
                            self.stdout.write(f"警告: 基因 {gene_id} 验证失败: {ve}")
                            
                            # 使用备用的symbol
                            gene_data['symbol'] = f"GENE_{gene_id}"
                            try:
                                obj, created = KeggGene.objects.get_or_create(
                                    gene_id=gene_id,
                                    organism=org_obj,
                                    defaults=gene_data
                                )
                                gene_objs[gene_id] = obj
                                imported_count += 1
                            except Exception:
                                continue
                        
                    except Exception as e:
                        self.stdout.write(f"错误: 导入基因 {gene.get('gene_id', 'N/A')} 失败: {e}")
                        continue
                
                self.stdout.write(f"成功导入 {imported_count} 个基因，处理了 {duplicate_count} 个重复符号")
                
                # 3. 导入关联（如果有关联数据）
                relation_file = self.data_dir / f"relations_{organism}.json"
                if relation_file.exists() and not self.skip_relations:
                    self.stdout.write("导入基因-通路关联...")
                    
                    with open(relation_file, 'r') as f:
                        relations = json.load(f)
                    
                    # 过滤有效的关联
                    valid_relations = []
                    for rel in relations:
                        if (rel.get('gene_id') in gene_objs and 
                            rel.get('pathway_id') in pathway_objs):
                            valid_relations.append(rel)
                    
                    self.stdout.write(f"有效关联数: {len(valid_relations)}/{len(relations)}")
                    
                    # 批量创建关联对象
                    pathway_gene_counts = {}
                    gene_pathway_counts = {}
                    
                    for rel in tqdm(valid_relations, desc="创建关联"):
                        gene_id = rel.get('gene_id')
                        pathway_id = rel.get('pathway_id')
                        
                        try:
                            # 检查关联是否已存在
                            if not KeggGenePathway.objects.filter(
                                gene=gene_objs[gene_id],
                                pathway=pathway_objs[pathway_id]
                            ).exists():
                                KeggGenePathway.objects.create(
                                    gene=gene_objs[gene_id],
                                    pathway=pathway_objs[pathway_id],
                                    relation_type='member',
                                    evidence='KEGG'
                                )
                            
                            # 计数
                            pathway_gene_counts[pathway_id] = \
                                pathway_gene_counts.get(pathway_id, 0) + 1
                            gene_pathway_counts[gene_id] = \
                                gene_pathway_counts.get(gene_id, 0) + 1
                        except Exception as e:
                            # 忽略关联创建错误，继续处理
                            continue
                    
                    # 4. 更新计数
                    self.stdout.write("更新计数...")
                    
                    # 更新通路基因计数
                    for pathway_id, count in pathway_gene_counts.items():
                        if pathway_id in pathway_objs:
                            pathway = pathway_objs[pathway_id]
                            pathway.gene_count = count
                            pathway.save()
                    
                    # 更新基因通路计数
                    for gene_id, count in gene_pathway_counts.items():
                        if gene_id in gene_objs:
                            gene = gene_objs[gene_id]
                            gene.pathway_count = count
                            gene.save()
                    
                    # 打印统计信息
                    self.stdout.write(self.style.SUCCESS(
                        f"\n导入统计:\n"
                        f"  生物体: {organism}\n"
                        f"  通路: {len(pathway_objs)}/{len(pathways)}\n"
                        f"  基因: {len(gene_objs)}/{len(genes)}\n"
                        f"  关联: {len(valid_relations)}/{len(relations)}\n"
                    ))
                else:
                    self.stdout.write("跳过关联导入或没有关联数据")
                    
        except Exception as e:
            self.stdout.write(self.style.ERROR(f"导入数据库失败: {e}"))
            import traceback
            traceback.print_exc()
            raise
    
    # 辅助方法
    def get_organism_name(self, code):
        """获取生物体名称"""
        organism_names = {
            'hsa': 'Homo sapiens (human)',
            'mmu': 'Mus musculus (mouse)',
            'rno': 'Rattus norvegicus (rat)',
            'dre': 'Danio rerio (zebrafish)',
            'cel': 'Caenorhabditis elegans',
            'dme': 'Drosophila melanogaster',
            'sce': 'Saccharomyces cerevisiae',
            'eco': 'Escherichia coli',
        }
        return organism_names.get(code, f"Organism {code}")
    
    def get_scientific_name(self, code):
        """获取学名"""
        scientific_names = {
            'hsa': 'Homo sapiens',
            'mmu': 'Mus musculus',
            'rno': 'Rattus norvegicus',
            'dre': 'Danio rerio',
            'cel': 'Caenorhabditis elegans',
            'dme': 'Drosophila melanogaster',
            'sce': 'Saccharomyces cerevisiae',
            'eco': 'Escherichia coli K-12 MG1655',
        }
        return scientific_names.get(code, '')
    
    def categorize_pathway(self, pathway_id, pathway_name):
        """根据通路ID和名称分类"""
        pathway_name_lower = pathway_name.lower()
        
        # 基于名称的分类
        if 'metabol' in pathway_name_lower:
            return 'metabolism'
        elif any(x in pathway_name_lower for x in ['signal', 'signaling', 'membrane', 'environment']):
            return 'environmental'
        elif any(x in pathway_name_lower for x in ['transcription', 'translation', 'replication', 'genetic']):
            return 'genetic'
        elif any(x in pathway_name_lower for x in ['cell', 'cycle', 'death', 'apoptosis', 'transport', 'cellular']):
            return 'cellular'
        elif any(x in pathway_name_lower for x in ['immune', 'endocrine', 'nervous', 'digestive', 'organismal']):
            return 'organismal'
        elif any(x in pathway_name_lower for x in ['cancer', 'disease', 'infection']):
            return 'human'
        elif any(x in pathway_name_lower for x in ['drug']):
            return 'drug'
        else:
            return 'genetic'