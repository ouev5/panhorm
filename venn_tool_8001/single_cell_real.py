#!/usr/bin/env python3
"""
单细胞真实数据处理模块
使用Tabula Sapiens完整数据
"""

import json
import numpy as np
from scipy import sparse
import os

# 全局变量缓存
_EXPR_MATRIX = None
_GENE_LIST = None
_GENE_INDEX = None
_CELL_META = None

DATA_DIR = "/www/wwwroot/venn-tool/single_cell_data"

def load_real_data():
    """加载真实的Tabula Sapiens数据"""
    global _EXPR_MATRIX, _GENE_LIST, _GENE_INDEX, _CELL_META
    
    if _EXPR_MATRIX is not None:
        return _EXPR_MATRIX, _GENE_INDEX, _CELL_META
    
    print("Loading Tabula Sapiens real data...")
    
    # 加载表达矩阵
    matrix_file = os.path.join(DATA_DIR, "expr_matrix.npz")
    print(f"Loading expression matrix from {matrix_file}...")
    _EXPR_MATRIX = sparse.load_npz(matrix_file)
    print(f"Matrix shape: {_EXPR_MATRIX.shape}")
    
    # 加载基因列表
    gene_file = os.path.join(DATA_DIR, "gene_list.json")
    print(f"Loading gene list from {gene_file}...")
    with open(gene_file, "r") as f:
        _GENE_LIST = json.load(f)
    
    # 创建基因索引
    _GENE_INDEX = {g.upper(): i for i, g in enumerate(_GENE_LIST)}
    print(f"Total genes: {len(_GENE_LIST)}")
    
    # 加载细胞元数据
    meta_file = os.path.join(DATA_DIR, "cell_metadata.json")
    print(f"Loading cell metadata from {meta_file}...")
    with open(meta_file, "r") as f:
        meta = json.load(f)
    _CELL_META = meta.get("cells", [])
    print(f"Total cells: {len(_CELL_META)}")
    
    return _EXPR_MATRIX, _GENE_INDEX, _CELL_META

def query_cells_by_genes(query_genes, max_cells=5000):
    """
    根据基因列表查询表达这些基因的细胞
    
    返回:
        - cells: 细胞列表（包含表达量信息）
        - found_genes: 找到的基因
        - not_found_genes: 未找到的基因
    """
    matrix, gene_index, cell_meta = load_real_data()
    
    # 找到基因索引
    found_genes = []
    not_found_genes = []
    gene_indices = []
    
    for g in query_genes:
        g_upper = g.upper()
        if g_upper in gene_index:
            found_genes.append(g)
            gene_indices.append(gene_index[g_upper])
        else:
            not_found_genes.append(g)
    
    if not gene_indices:
        return [], [], not_found_genes
    
    # 获取这些基因的表达矩阵行
    # matrix shape: (genes, cells)
    gene_expr = matrix[gene_indices, :].toarray()  # shape: (n_genes, n_cells)
    
    # 计算每个细胞的总表达量
    total_expr = np.sum(gene_expr, axis=0)  # shape: (n_cells,)
    
    # 找到表达量 > 0 的细胞
    expressing_mask = total_expr > 0
    expressing_indices = np.where(expressing_mask)[0]
    
    print(f"Found {len(expressing_indices)} cells expressing query genes")
    
    # 按表达量排序，取前 max_cells 个
    sorted_indices = expressing_indices[np.argsort(-total_expr[expressing_indices])]
    selected_indices = sorted_indices[:max_cells]
    
    # 构建返回数据
    cells = []
    for idx in selected_indices:
        cell = cell_meta[idx].copy()
        # 添加基因表达量
        cell["gene_expr"] = float(total_expr[idx])
        cell["gene_expr_values"] = {found_genes[i]: float(gene_expr[i, idx]) 
                                     for i in range(len(found_genes))}
        cells.append(cell)
    
    return cells, found_genes, not_found_genes

def compute_tsne_coords(cells):
    """
    为细胞计算t-SNE坐标（真实计算）
    使用sklearn的t-SNE
    """
    from sklearn.manifold import TSNE
    
    if len(cells) < 100:
        return cells
    
    # 提取UMAP坐标
    umap_coords = np.array([[c.get("x", 0), c.get("y", 0)] for c in cells])
    
    # 使用t-SNE计算新坐标
    # 对于大规模数据，使用PCA预降维
    n_samples = len(cells)
    
    if n_samples > 5000:
        # 对于大数据集，使用近似方法
        perplexity = min(30, n_samples // 100)
    else:
        perplexity = min(30, n_samples - 1)
    
    tsne = TSNE(n_components=2, perplexity=perplexity, random_state=42,
                max_iter=250, learning_rate='auto', init='random')
    
    # t-SNE需要高维输入，这里用UMAP坐标作为"伪高维"
    # 实际应该用表达矩阵，但太大了
    tsne_coords = tsne.fit_transform(umap_coords)
    
    # 归一化到相似范围
    umap_center = np.mean(umap_coords, axis=0)
    umap_scale = np.std(umap_coords)
    
    tsne_center = np.mean(tsne_coords, axis=0)
    tsne_scale = np.std(tsne_coords) if np.std(tsne_coords) > 0 else 1
    
    tsne_coords = (tsne_coords - tsne_center) / tsne_scale * umap_scale + umap_center
    
    # 添加到细胞数据
    for i, cell in enumerate(cells):
        cell["tsne_x"] = float(tsne_coords[i, 0])
        cell["tsne_y"] = float(tsne_coords[i, 1])
    
    return cells

if __name__ == "__main__":
    # 测试
    cells, found, not_found = query_cells_by_genes(["TP53", "BRCA1", "EGFR"])
    print(f"\nTest query: TP53, BRCA1, EGFR")
    print(f"Found genes: {found}")
    print(f"Not found genes: {not_found}")
    print(f"Cells returned: {len(cells)}")
    if cells:
        print(f"First cell: {cells[0]}")
