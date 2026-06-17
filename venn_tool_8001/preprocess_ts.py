import anndata as ad
import json
import os

print("=== 加载Tabula Sapiens数据 ===")
adata = ad.read_h5ad("single_cell_data/tabula_sapiens_chunk0.h5ad")
print(f"细胞数: {adata.n_obs}")
print(f"基因数: {adata.n_vars}")

# 创建基因名到Ensembl ID的映射
gene_name_to_ensg = {}
for ensg in adata.var_names:
    gene_name = adata.var.loc[ensg, 'gene_name']
    if gene_name and isinstance(gene_name, str):
        gene_name_to_ensg[gene_name.upper()] = ensg

print(f"\n基因映射数: {len(gene_name_to_ensg)}")

# 提取UMAP坐标和元数据
import numpy as np

# UMAP坐标
umap_coords = adata.obsm['X_umap']
print(f"UMAP形状: {umap_coords.shape}")

# 细胞类型
cell_types = adata.obs['cell_ontology_class'].fillna('Unknown').tolist()
unique_types = list(set(cell_types))
print(f"细胞类型数: {len(unique_types)}")

# 组织
tissues = adata.obs['organ_tissue'].fillna('Unknown').tolist()
unique_tissues = list(set(tissues))
print(f"组织类型数: {len(unique_tissues)}")

# 创建颜色映射
import colorsys
def get_colors(n):
    colors = []
    for i in range(n):
        hue = i / n
        rgb = colorsys.hsv_to_rgb(hue, 0.7, 0.8)
        color = '#{:02x}{:02x}{:02x}'.format(int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
        colors.append(color)
    return colors

type_colors = get_colors(len(unique_types))
type_color_map = {t: type_colors[i] for i, t in enumerate(unique_types)}

tissue_colors = get_colors(len(unique_tissues))
tissue_color_map = {t: tissue_colors[i] for i, t in enumerate(unique_tissues)}

# 保存预处理数据
print("\n保存预处理数据...")

# 保存UMAP坐标和元数据（用于前端可视化）
cell_data = []
for i in range(adata.n_obs):
    cell_data.append({
        'x': float(umap_coords[i, 0]),
        'y': float(umap_coords[i, 1]),
        'cell_type': cell_types[i],
        'tissue': tissues[i],
        'type_color': type_color_map[cell_types[i]],
        'tissue_color': tissue_color_map[tissues[i]]
    })

# 保存为JSON
with open('single_cell_data/cell_metadata.json', 'w') as f:
    json.dump({
        'cells': cell_data[:50000],  # 限制50000个细胞避免文件过大
        'cell_types': unique_types,
        'tissues': unique_tissues,
        'type_colors': type_color_map,
        'tissue_colors': tissue_color_map
    }, f)

print(f"保存细胞元数据: {len(cell_data)} 个细胞")

# 保存基因映射
with open('single_cell_data/gene_mapping.json', 'w') as f:
    json.dump(gene_name_to_ensg, f)

print(f"保存基因映射: {len(gene_name_to_ensg)} 个基因")

# 保存表达矩阵的稀疏表示（用于查询特定基因的表达）
from scipy import sparse
if not sparse.issparse(adata.X):
    adata.X = sparse.csr_matrix(adata.X)

# 创建基因索引
gene_to_idx = {gene_name_to_ensg[k]: i for i, k in enumerate(gene_name_to_ensg.keys()) if gene_name_to_ensg[k] in adata.var_names}
print(f"\n基因索引数: {len(gene_to_idx)}")

# 保存adata对象供后端API使用
print("\n保存处理后的AnnData...")
adata.write_h5ad('single_cell_data/ts_processed.h5ad')
print("完成!")

print(f"\n=== 最终文件 ===")
for f in ['single_cell_data/cell_metadata.json', 'single_cell_data/gene_mapping.json', 'single_cell_data/ts_processed.h5ad']:
    if os.path.exists(f):
        size_mb = os.path.getsize(f) / 1024 / 1024
        print(f"  {f}: {size_mb:.1f} MB")
