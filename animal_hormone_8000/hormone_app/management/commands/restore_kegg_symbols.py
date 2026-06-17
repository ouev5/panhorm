# management/commands/restore_kegg_symbols.py
import json
from pathlib import Path
from django.core.management.base import BaseCommand
from django.db import transaction
from hormone_app.models import KeggGene, KeggOrganism
from django.core.exceptions import ValidationError

class Command(BaseCommand):
    help = '从KEGG原始数据恢复正确的基因符号'
    
    def add_arguments(self, parser):
        parser.add_argument('--organism', type=str, default='hsa', help='生物体代码')
        parser.add_argument('--data-dir', type=str, default='data/kegg', help='数据目录')
        parser.add_argument('--dry-run', action='store_true', help='模拟运行')
    
    def handle(self, *args, **options):
        organism = options['organism']
        data_dir = Path(options['data_dir'])
        dry_run = options['dry_run']
        
        try:
            org = KeggOrganism.objects.get(code=organism)
            
            # 读取原始基因数据文件
            gene_file = data_dir / f"genes_{organism}.json"
            if not gene_file.exists():
                self.stdout.write(self.style.ERROR(f"基因数据文件不存在: {gene_file}"))
                return
            
            with open(gene_file, 'r') as f:
                raw_genes = json.load(f)
            
            self.stdout.write(f"从原始数据加载了 {len(raw_genes)} 个基因记录")
            
            # 先统计所有可能重复的符号
            symbol_counter = {}
            for raw_gene in raw_genes:
                gene_id = raw_gene['gene_id']
                if 'symbol' in raw_gene and raw_gene['symbol'] and not raw_gene['symbol'].startswith('GENE_'):
                    symbol = raw_gene['symbol']
                    symbol_counter[symbol] = symbol_counter.get(symbol, 0) + 1
                elif 'name' in raw_gene:
                    name_parts = raw_gene['name'].split(';')
                    if name_parts:
                        symbol = name_parts[0].strip()
                        if symbol and not symbol.startswith('GENE_'):
                            symbol_counter[symbol] = symbol_counter.get(symbol, 0) + 1
            
            self.stdout.write(f"发现 {len(symbol_counter)} 个唯一符号")
            
            # 找出重复的符号
            duplicate_symbols = {s: c for s, c in symbol_counter.items() if c > 1}
            self.stdout.write(f"重复符号数: {len(duplicate_symbols)}")
            
            # 为每个重复符号创建计数
            symbol_usage = {}
            
            with transaction.atomic():
                updated = 0
                skipped = 0
                error_count = 0
                
                for raw_gene in raw_genes:
                    try:
                        gene_id = raw_gene['gene_id']
                        
                        # 获取正确的符号
                        correct_symbol = None
                        if 'symbol' in raw_gene and raw_gene['symbol'] and not raw_gene['symbol'].startswith('GENE_'):
                            correct_symbol = raw_gene['symbol']
                        elif 'name' in raw_gene:
                            name_parts = raw_gene['name'].split(';')
                            if name_parts:
                                symbol = name_parts[0].strip()
                                if symbol and not symbol.startswith('GENE_'):
                                    correct_symbol = symbol
                        
                        if not correct_symbol:
                            continue
                        
                        # 处理重复符号
                        if correct_symbol in duplicate_symbols:
                            # 增加计数器
                            symbol_usage[correct_symbol] = symbol_usage.get(correct_symbol, 0) + 1
                            
                            # 如果不是第一个使用此符号的基因，添加后缀
                            if symbol_usage[correct_symbol] > 1:
                                correct_symbol = f"{correct_symbol}_{symbol_usage[correct_symbol]}"
                        
                        # 限制长度
                        if len(correct_symbol) > 200:
                            correct_symbol = correct_symbol[:197] + "..."
                        
                        try:
                            gene = KeggGene.objects.get(gene_id=gene_id, organism=org)
                            
                            if gene.symbol != correct_symbol:
                                if not dry_run:
                                    # 直接更新，跳过模型验证
                                    KeggGene.objects.filter(
                                        gene_id=gene_id, 
                                        organism=org
                                    ).update(symbol=correct_symbol)
                                
                                updated += 1
                                if updated <= 50:
                                    self.stdout.write(f"  ✓ {gene_id}: {gene.symbol[:30]} -> {correct_symbol[:30]}")
                            else:
                                skipped += 1
                                
                        except KeggGene.DoesNotExist:
                            continue
                            
                    except Exception as e:
                        error_count += 1
                        if error_count <= 10:
                            self.stdout.write(self.style.WARNING(f"  错误 {gene_id}: {e}"))
                        continue
                
                self.stdout.write(f"\n统计:")
                self.stdout.write(f"  - 更新: {updated}")
                self.stdout.write(f"  - 跳过: {skipped}")
                self.stdout.write(f"  - 错误: {error_count}")
                
                if dry_run:
                    self.stdout.write(self.style.WARNING("这是模拟运行，未实际修改数据库"))
                else:
                    self.stdout.write(self.style.SUCCESS(f"成功更新 {updated} 个基因符号"))
            
        except KeggOrganism.DoesNotExist:
            self.stdout.write(self.style.ERROR(f"生物体 {organism} 不存在"))
        except Exception as e:
            self.stdout.write(self.style.ERROR(f"错误: {e}"))
            import traceback
            traceback.print_exc()