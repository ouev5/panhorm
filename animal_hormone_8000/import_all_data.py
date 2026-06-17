import os
import django
import pandas as pd
from django.conf import settings

# 配置 Django 环境（必须）
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "animal_hormone.settings")  # 替换为你的项目settings路径
django.setup()

# 导入需要操作的模型（替换为你的应用名）
from hormone_app.models import (
    HormoneRelatedGene,
    HormoneReceptorInfo,
    HormoneReceptorFull
)


def import_hormone_related_gene(file_path):
    """导入 HormoneRelatedGene (对应：激素-基因-疾病-DO1.xlsx)"""
    try:
        # 读取 Excel/CSV 文件（自动识别格式）
        if file_path.endswith('.xlsx'):
            df = pd.read_excel(file_path, dtype=str, keep_default_na=False)
        else:
            df = pd.read_csv(file_path, dtype=str, keep_default_na=False)

        # 遍历数据，批量导入
        batch = []
        for _, row in df.iterrows():
            # 处理空值（将"NA"或空字符串转为None）
            def get_val(key):
                val = row.get(key, "")
                return val if val not in ("NA", "") else None

            # 构造模型对象
            obj = HormoneRelatedGene(
                id=get_val('ID'),  # 注意：原数据ID是主键，需确保唯一
                hormone_name=get_val('hormone_name'),
                related_genes=get_val('related_genes'),
                pmid=int(get_val('pmid')) if get_val('pmid') else None,
                gene_sequence=get_val('gene_sequence'),
                related_diseases=get_val('related_diseases'),
                DO_Match_disease_names=get_val('DO_Match_disease_names'),
                DO_ID=get_val('DO_ID'),
                DO_Standardized_Terminology=get_val('DO_Standardized_Terminology')
            )
            batch.append(obj)

            # 每1000条批量插入（提高效率）
            if len(batch) >= 1000:
                HormoneRelatedGene.objects.bulk_create(batch, ignore_conflicts=True)  # 忽略主键冲突
                batch = []

        # 插入剩余数据
        if batch:
            HormoneRelatedGene.objects.bulk_create(batch, ignore_conflicts=True)
        print(f"✅ HormoneRelatedGene 导入完成，共导入 {len(df)} 条数据")

    except Exception as e:
        print(f"❌ HormoneRelatedGene 导入失败：{str(e)}")


def import_hormone_receptor_info(file_path):
    """导入 HormoneReceptorInfo (对应：非肽类激素受体_dedup.csv)"""
    try:
        df = pd.read_csv(file_path, dtype=str, keep_default_na=False)
        batch = []
        for _, row in df.iterrows():
            def get_val(key):
                val = row.get(key, "")
                return val if val not in ("NA", "") else None

            obj = HormoneReceptorInfo(
                id=get_val('ID'),
                pubchem_id=get_val('Pubchem ID'),
                hormone_name=get_val('hormone name'),
                receptor_uniprot_id=get_val('receptor uniProt ID'),
                receptor_coding_genes=get_val('receptor coding genes'),
                receptor_name=get_val('receptor name'),
                receptor_species_name=get_val('receptor species name'),
                receptor_coding_genes_sequence=get_val('receptor coding genes sequence')
            )
            batch.append(obj)
            if len(batch) >= 1000:
                HormoneReceptorInfo.objects.bulk_create(batch, ignore_conflicts=True)
                batch = []
        if batch:
            HormoneReceptorInfo.objects.bulk_create(batch, ignore_conflicts=True)
        print(f"✅ HormoneReceptorInfo 导入完成，共导入 {len(df)} 条数据")

    except Exception as e:
        print(f"❌ HormoneReceptorInfo 导入失败：{str(e)}")


def import_hormone_receptor_full(file_path):
    """导入 HormoneReceptorFull (对应：肽类激素受体_dedup.csv)"""
    try:
        df = pd.read_csv(file_path, dtype=str, keep_default_na=False)
        batch = []
        for _, row in df.iterrows():
            def get_val(key):
                val = row.get(key, "")
                return val if val not in ("NA", "") else None

            obj = HormoneReceptorFull(
                id=get_val('ID'),
                hormone_uniprot_id=get_val('hormone uniProt ID'),
                hormone_name=get_val('hormone name'),
                hormone_species_name=get_val('hormone species name'),
                hormone_coding_genes=get_val('hormone coding genes'),
                hormone_coding_genes_sequence=get_val('hormone coding genes sequence'),
                receptor_uniprot_id=get_val('receptor uniProt ID'),
                receptor_name=get_val('receptor name'),
                receptor_species_name=get_val('receptor species name'),
                receptor_coding_genes=get_val('receptor coding genes'),
                receptor_coding_genes_sequence=get_val('receptor coding genes sequence')
            )
            batch.append(obj)
            if len(batch) >= 1000:
                HormoneReceptorFull.objects.bulk_create(batch, ignore_conflicts=True)
                batch = []
        if batch:
            HormoneReceptorFull.objects.bulk_create(batch, ignore_conflicts=True)
        print(f"✅ HormoneReceptorFull 导入完成，共导入 {len(df)} 条数据")

    except Exception as e:
        print(f"❌ HormoneReceptorFull 导入失败：{str(e)}")


if __name__ == "__main__":
    # 你的文件路径（对应截图中的三个文件）
    file_paths = {
        "HormoneRelatedGene": "/www/wwwroot/default/animal_hormone/media/csv/激素-基因-疾病-DO1.xlsx",
        "HormoneReceptorInfo": "/www/wwwroot/default/animal_hormone/media/csv/非肽类激素受体_dedup.csv",
        "HormoneReceptorFull": "/www/wwwroot/default/animal_hormone/media/csv/肽类激素受体_dedup.csv"
    }

    # 依次导入三个表
    import_hormone_related_gene(file_paths["HormoneRelatedGene"])
    import_hormone_receptor_info(file_paths["HormoneReceptorInfo"])
    import_hormone_receptor_full(file_paths["HormoneReceptorFull"])