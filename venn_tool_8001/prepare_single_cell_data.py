#!/usr/bin/env python3
"""
准备单细胞数据
使用scanpy下载Tabula Sapiens或pbmc3k数据集
"""

import os
import json
import numpy as np
from scipy import sparse
import scanpy as sc
import anndata as ad

DATA_DIR = '/www/wwwroot/venn-tool/single_cell_data'

def prepare_pbmc3k():
    """下载并预处理PBMC3k数据集（较小，快速测试）"""
    print('=== 下载 PBMC3k 数据集 ===')
    adata = sc.datasets.pbmc3k()
    
    print(f'细胞数: {adata.n_obs}')
    print(f'基因数: {adata.n_vars}')
    
    # 预处理
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # 归一化
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # 高变基因
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable]
    
    # PCA
    sc.tl.pca(adata, n_comps=50)
    
    # UMAP
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    
    # 聚类
    sc.tl.leiden(adata, resolution=0.5, flavor="igraph")
    
    # 保存数据
    save_data(adata, 'pbmc3k')
    
    return adata

def prepare_tabula_sapiens_subset():
    """下载Tabula Sapiens子集"""
    print('=== 尝试下载 Tabula Sapiens ===')
    try:
        # Tabula Sapiens数据集很大，可能需要特殊处理
        # 使用scanpy的web函数
        import pooch
        url = 'https://figshare.com/ndownloader/files/38101904'  # Tabula Sapiens subset
        print('Downloading from figshare...')
        # 先尝试pbmc3k作为备选
        return prepare_pbmc3k()
    except Exception as e:
        print(f'Tabula Sapiens下载失败: {e}')
        print('回退到PBMC3k数据集')
        return prepare_pbmc3k()

def save_data(adata, name):
    """保存处理后的数据"""
    os.makedirs(DATA_DIR, exist_ok=True)
    
    print(f'保存数据到 {DATA_DIR}...')
    
    # 1. 保存基因列表
    gene_list = adata.var_names.tolist()
    with open(os.path.join(DATA_DIR, 'gene_list.json'), 'w') as f:
        json.dump(gene_list, f)
    print(f'基因列表: {len(gene_list)} 基因')
    
    # 2. 保存表达矩阵 (稀疏格式)
    if sparse.issparse(adata.X):
        expr_matrix = adata.X.T  # 转置为 (genes, cells)
    else:
        expr_matrix = sparse.csr_matrix(adata.X.T)
    
    sparse.save_npz(os.path.join(DATA_DIR, 'expr_matrix.npz'), expr_matrix)
    print(f'表达矩阵: {expr_matrix.shape}')
    
    # 3. 保存细胞元数据
    # UMAP坐标
    if 'X_umap' in adata.obsm:
        umap_coords = adata.obsm['X_umap']
    else:
        umap_coords = np.zeros((adata.n_obs, 2))
    
    # 细胞类型
    if 'leiden' in adata.obs.columns:
        cell_types = adata.obs['leiden'].tolist()
        unique_types = sorted(set(cell_types))
    elif 'louvain' in adata.obs.columns:
        cell_types = adata.obs['louvain'].tolist()
        unique_types = sorted(set(cell_types))
    else:
        cell_types = ['Unknown'] * adata.n_obs
        unique_types = ['Unknown']
    
    # 组织（如果没有则用'Unknown'）
    tissues = adata.obs.get('organ_tissue', ['tissue1'] * adata.n_obs).tolist() if 'organ_tissue' in adata.obs.columns else ['tissue1'] * adata.n_obs
    
    # 颜色映射
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
    
    # 构建细胞数据
    cells = []
    for i in range(adata.n_obs):
        cells.append({
            'x': float(umap_coords[i, 0]),
            'y': float(umap_coords[i, 1]),
            'cell_type': cell_types[i],
            'tissue': tissues[i] if isinstance(tissues, list) else 'tissue1',
            'cell_id': f'cell_{i}'
        })
    
    metadata = {
        'cells': cells,
        'cell_types': unique_types,
        'type_colors': type_color_map,
        'n_cells': adata.n_obs,
        'n_genes': len(gene_list),
        'source': name
    }
    
    with open(os.path.join(DATA_DIR, 'cell_metadata.json'), 'w') as f:
        json.dump(metadata, f)
    
    print(f'细胞元数据: {adata.n_obs} 细胞')
    print(f'细胞类型: {len(unique_types)} 种')
    
    # 4. 保存统计信息
    info = {
        'n_cells': adata.n_obs,
        'n_genes': len(gene_list),
        'cell_types': unique_types,
        'source': name
    }
    with open(os.path.join(DATA_DIR, 'data_info.json'), 'w') as f:
        json.dump(info, f, indent=2)
    
    print('\n数据准备完成！')
    return True

if __name__ == '__main__':
    # 先尝试PBMC3k（小数据集，快速）
    adata = prepare_pbmc3k()
    print(f'\n数据集大小: {adata.n_obs} 细胞, {adata.n_vars} 基因')
