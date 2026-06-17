import csv
from django.core.management.base import BaseCommand
from hormone_app.models import (
    HormoneRelatedGene,
    HormoneReceptorInfo,
    HormoneReceptorFull
)

class Command(BaseCommand):
    help = '导入三个文件的数据到数据库，保留NA值作为字符串，已存在则删除更新'

    def add_arguments(self, parser):
        parser.add_argument('file1', type=str, help='第一个文件路径（hormone_related_genes.csv）')
        parser.add_argument('file2', type=str, help='第二个文件路径（hormone_receptor_info.csv）')
        parser.add_argument('file3', type=str, help='第三个文件路径（hormone_receptor_full.csv）')

    def handle(self, *args, **options):
        try:
            # 导入三个文件并统计导入数量
            count1 = self.import_file1(options['file1'])
            count2 = self.import_file2(options['file2'])
            count3 = self.import_file3(options['file3'])
            
            self.stdout.write(self.style.SUCCESS(
                f'数据导入完成！\n'
                f'文件1导入: {count1} 条记录\n'
                f'文件2导入: {count2} 条记录\n'
                f'文件3导入: {count3} 条记录\n'
                f'(NA值已保留，已存在数据已更新)'
            ))
        except Exception as e:
            self.stderr.write(self.style.ERROR(f'导入过程出错: {str(e)}'))

    def preserve_na(self, value):
        """统一处理NA值，空字符串或仅空格的情况转为'NA'，保留原始NA"""
        if value is None:
            return 'NA'
        # 去除首尾空格
        stripped = value.strip()
        # 空值或仅空格的情况
        if not stripped:
            return 'NA'
        # 已经是NA的情况保持不变
        return stripped

    def safe_int(self, value, field_name, row_num):
        """安全转换为整数，处理空值和NA，提供更详细的错误信息"""
        if value is None:
            self.stderr.write(self.style.WARNING(f"第 {row_num} 行 '{field_name}' 为空，跳过该行"))
            return None
            
        stripped = value.strip()
        if not stripped:
            self.stderr.write(self.style.WARNING(f"第 {row_num} 行 '{field_name}' 为空，跳过该行"))
            return None
            
        if stripped.upper() == 'NA':
            self.stderr.write(self.style.WARNING(f"第 {row_num} 行 '{field_name}' 为NA，无法转换为整数，跳过该行"))
            return None
            
        try:
            return int(stripped)
        except ValueError:
            self.stderr.write(self.style.WARNING(f"第 {row_num} 行 '{field_name}' 值 '{stripped}' 不是有效整数，跳过该行"))
            return None

    def import_file1(self, file_path):
        """导入hormone_related_genes.csv数据"""
        imported = 0
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                # 验证CSV列是否完整
                required_columns = ['ID', 'hormone_name', 'related_genes', 'pmid', 'gene_sequence']
                for col in required_columns:
                    if col not in reader.fieldnames:
                        raise ValueError(f"文件 {file_path} 缺少必要列: {col}")

                for idx, row in enumerate(reader, start=1):
                    # 跳过全空行
                    if all(not value.strip() for value in row.values()):
                        continue

                    # 处理ID
                    row_id = self.safe_int(row.get('ID'), 'ID', idx)
                    if row_id is None:
                        continue

                    # 处理pmid
                    pmid = self.safe_int(row.get('pmid'), 'pmid', idx)
                    if pmid is None:
                        continue

                    # 检查数据是否已存在，存在则删除
                    existing = HormoneRelatedGene.objects.filter(id=row_id).first()
                    if existing:
                        existing.delete()
                        self.stdout.write(self.style.WARNING(f"第 {idx} 行 ID {row_id} 已存在，已更新"))

                    # 处理其他字段，保留NA
                    hormone_name = self.preserve_na(row.get('hormone_name'))
                    related_genes = self.preserve_na(row.get('related_genes'))
                    gene_sequence = self.preserve_na(row.get('gene_sequence'))

                    # 创建记录
                    HormoneRelatedGene.objects.create(
                        id=row_id,
                        hormone_name=hormone_name,
                        related_genes=related_genes,
                        pmid=pmid,
                        gene_sequence=gene_sequence
                    )
                    imported += 1

                    # 每100行显示一次进度
                    if imported % 100 == 0:
                        self.stdout.write(f"文件1已导入 {imported} 条记录...")

            self.stdout.write(self.style.SUCCESS(f"文件1导入完成，共导入 {imported} 条记录"))
            return imported
            
        except FileNotFoundError:
            raise ValueError(f"文件 {file_path} 不存在")
        except Exception as e:
            raise ValueError(f"导入文件1时出错: {str(e)}")

    def import_file2(self, file_path):
        """导入hormone_receptor_info.csv数据"""
        imported = 0
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                required_columns = ['ID', 'hormone name', 'receptor uniProt ID', 
                                   'receptor coding genes', 'receptor name',
                                   'receptor species name', 'receptor coding genes sequence',
                                   'Pubchem ID']
                for col in required_columns:
                    if col not in reader.fieldnames:
                        raise ValueError(f"文件 {file_path} 缺少必要列: {col}")

                for idx, row in enumerate(reader, start=1):
                    # 跳过全空行
                    if all(not value.strip() for value in row.values()):
                        continue

                    row_id = self.safe_int(row.get('ID'), 'ID', idx)
                    if row_id is None:
                        continue

                    pubchem_id = self.safe_int(row.get('Pubchem ID'), 'Pubchem ID', idx)
                    if pubchem_id is None:
                        continue

                    # 检查数据是否已存在
                    existing = HormoneReceptorInfo.objects.filter(id=row_id).first()
                    if existing:
                        existing.delete()
                        self.stdout.write(self.style.WARNING(f"第 {idx} 行 ID {row_id} 已存在，已更新"))

                    # 处理其他字段
                    hormone_name = self.preserve_na(row.get('hormone name'))
                    receptor_uniprot_id = self.preserve_na(row.get('receptor uniProt ID'))
                    receptor_coding_genes = self.preserve_na(row.get('receptor coding genes'))
                    receptor_name = self.preserve_na(row.get('receptor name'))
                    receptor_species_name = self.preserve_na(row.get('receptor species name'))
                    receptor_coding_genes_sequence = self.preserve_na(row.get('receptor coding genes sequence'))

                    HormoneReceptorInfo.objects.create(
                        id=row_id,
                        pubchem_id=pubchem_id,
                        hormone_name=hormone_name,
                        receptor_uniprot_id=receptor_uniprot_id,
                        receptor_coding_genes=receptor_coding_genes,
                        receptor_name=receptor_name,
                        receptor_species_name=receptor_species_name,
                        receptor_coding_genes_sequence=receptor_coding_genes_sequence
                    )
                    imported += 1

                    if imported % 100 == 0:
                        self.stdout.write(f"文件2已导入 {imported} 条记录...")

            self.stdout.write(self.style.SUCCESS(f"文件2导入完成，共导入 {imported} 条记录"))
            return imported
            
        except FileNotFoundError:
            raise ValueError(f"文件 {file_path} 不存在")
        except Exception as e:
            raise ValueError(f"导入文件2时出错: {str(e)}")

    def import_file3(self, file_path):
        """导入hormone_receptor_full.csv数据"""
        imported = 0
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                required_columns = ['ID', 'hormone uniProt ID', 'hormone name',
                                   'hormone species name', 'hormone coding genes',
                                   'hormone coding genes sequence', 'receptor uniProt ID',
                                   'receptor name', 'receptor species name',
                                   'receptor coding genes', 'receptor coding genes sequence']
                for col in required_columns:
                    if col not in reader.fieldnames:
                        raise ValueError(f"文件 {file_path} 缺少必要列: {col}")

                for idx, row in enumerate(reader, start=1):
                    # 跳过全空行
                    if all(not value.strip() for value in row.values()):
                        continue

                    row_id = self.safe_int(row.get('ID'), 'ID', idx)
                    if row_id is None:
                        continue

                    # 检查数据是否已存在
                    existing = HormoneReceptorFull.objects.filter(id=row_id).first()
                    if existing:
                        existing.delete()
                        self.stdout.write(self.style.WARNING(f"第 {idx} 行 ID {row_id} 已存在，已更新"))

                    # 创建记录
                    HormoneReceptorFull.objects.create(
                        id=row_id,
                        hormone_uniprot_id=self.preserve_na(row.get('hormone uniProt ID')),
                        hormone_name=self.preserve_na(row.get('hormone name')),
                        hormone_species_name=self.preserve_na(row.get('hormone species name')),
                        hormone_coding_genes=self.preserve_na(row.get('hormone coding genes')),
                        hormone_coding_genes_sequence=self.preserve_na(row.get('hormone coding genes sequence')),
                        receptor_uniprot_id=self.preserve_na(row.get('receptor uniProt ID')),
                        receptor_name=self.preserve_na(row.get('receptor name')),
                        receptor_species_name=self.preserve_na(row.get('receptor species name')),
                        receptor_coding_genes=self.preserve_na(row.get('receptor coding genes')),
                        receptor_coding_genes_sequence=self.preserve_na(row.get('receptor coding genes sequence'))
                    )
                    imported += 1

                    if imported % 100 == 0:
                        self.stdout.write(f"文件3已导入 {imported} 条记录...")

            self.stdout.write(self.style.SUCCESS(f"文件3导入完成，共导入 {imported} 条记录"))
            return imported
            
        except FileNotFoundError:
            raise ValueError(f"文件 {file_path} 不存在")
        except Exception as e:
            raise ValueError(f"导入文件3时出错: {str(e)}")
