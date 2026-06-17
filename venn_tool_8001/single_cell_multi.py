#!/usr/bin/env python3
"""
通用单细胞数据处理模块 - 支持多物种
通过数据目录切换物种: single_cell_data (human) / single_cell_data_mouse (mouse)
"""
import json
import numpy as np
from scipy import sparse
import os

# 每个物种的缓存
_CACHES = {}

DATA_DIRS = {
    "human": "/www/wwwroot/venn-tool/single_cell_data",
    "mouse": "/www/wwwroot/venn-tool/single_cell_data_mouse",
}

def load_data(species="human"):
    """加载指定物种的单细胞数据"""
    if species in _CACHES:
        return _CACHES[species]
    
    data_dir = DATA_DIRS.get(species)
    if not data_dir or not os.path.exists(os.path.join(data_dir, "expr_matrix.npz")):
        return None, None, None
    
    print(f"Loading {species} single cell data from {data_dir}...")
    
    # 表达矩阵
    matrix_file = os.path.join(data_dir, "expr_matrix.npz")
    expr_matrix = sparse.load_npz(matrix_file)
    print(f"  Matrix: {expr_matrix.shape}")
    
    # 基因列表
    gene_file = os.path.join(data_dir, "gene_list.json")
    with open(gene_file, "r") as f:
        gene_list = json.load(f)
    gene_index = {g.upper(): i for i, g in enumerate(gene_list)}
    print(f"  Genes: {len(gene_list)}")
    
    # 细胞元数据
    meta_file = os.path.join(data_dir, "cell_metadata.json")
    with open(meta_file, "r") as f:
        meta = json.load(f)
    cell_meta = meta.get("cells", [])
    print(f"  Cells: {len(cell_meta)}")
    
    _CACHES[species] = (expr_matrix, gene_index, cell_meta)
    return expr_matrix, gene_index, cell_meta


def query_cells_by_genes(query_genes, species="human", max_cells=5000):
    """根据基因列表查询表达这些基因的细胞"""
    matrix, gene_index, cell_meta = load_data(species)
    
    if matrix is None:
        return [], [], query_genes
    
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
        return [], found_genes, not_found_genes
    
    gene_expr = matrix[gene_indices, :].toarray()
    total_expr = np.sum(gene_expr, axis=0)
    
    expressing_mask = total_expr > 0
    expressing_indices = np.where(expressing_mask)[0]
    
    print(f"Found {len(expressing_indices)} cells expressing query genes ({species})")
    
    sorted_indices = expressing_indices[np.argsort(-total_expr[expressing_indices])]
    selected_indices = sorted_indices[:max_cells]
    
    cells = []
    for idx in selected_indices:
        cell = cell_meta[idx].copy()
        cell["gene_expr"] = float(total_expr[idx])
        cell["gene_expr_values"] = {found_genes[i]: float(gene_expr[i, idx]) 
                                     for i in range(len(found_genes))}
        cells.append(cell)
    
    return cells, found_genes, not_found_genes


def compute_tsne_coords(cells):
    """为细胞计算t-SNE坐标"""
    from sklearn.manifold import TSNE
    
    if len(cells) < 100:
        return cells
    
    umap_coords = np.array([[c.get("x", 0), c.get("y", 0)] for c in cells])
    n_samples = len(cells)
    perplexity = min(30, (n_samples - 1) // 1)
    if n_samples > 5000:
        perplexity = min(30, n_samples // 100)
    else:
        perplexity = min(30, n_samples - 1)
    
    tsne = TSNE(n_components=2, perplexity=perplexity, random_state=42,
                max_iter=250, learning_rate='auto', init='random')
    tsne_coords = tsne.fit_transform(umap_coords)
    
    umap_center = np.mean(umap_coords, axis=0)
    umap_scale = np.std(umap_coords)
    tsne_center = np.mean(tsne_coords, axis=0)
    tsne_scale = np.std(tsne_coords) if np.std(tsne_coords) > 0 else 1
    tsne_coords = (tsne_coords - tsne_center) / tsne_scale * umap_scale + umap_center
    
    for i, cell in enumerate(cells):
        cell["tsne_x"] = float(tsne_coords[i, 0])
        cell["tsne_y"] = float(tsne_coords[i, 1])
    
    return cells


if __name__ == "__main__":
    print("=== Testing Human ===")
    cells, found, nf = query_cells_by_genes(["TP53", "BRCA1", "EGFR"], species="human")
    print(f"Found: {found}, Not found: {nf}, Cells: {len(cells)}")
    
    print("\n=== Testing Mouse ===")
    cells, found, nf = query_cells_by_genes(["Trp53", "Brca1", "Egfr"], species="mouse")
    print(f"Found: {found}, Not found: {nf}, Cells: {len(cells)}")
    if cells:
        print(f"First: {cells[0]['cell_type']} @ {cells[0]['tissue']}")
