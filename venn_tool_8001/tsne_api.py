#!/usr/bin/env python3
"""
t-SNE 计算模块
支持实时计算和预计算两种模式
"""

import numpy as np
from sklearn.manifold import TSNE
import os
import json
import logging

logger = logging.getLogger(__name__)

# 全局缓存
_tsne_cache = {
    'umap_coords': None,
    'tsne_coords': None,
    'metadata': None,
    'gene_list': None
}

def load_single_cell_data(base_path='/www/wwwroot/venn-tool/single_cell_data'):
    """加载单细胞数据"""
    metadata_path = os.path.join(base_path, 'cell_metadata.json')
    gene_list_path = os.path.join(base_path, 'gene_list.json')
    expr_path = os.path.join(base_path, 'expr_matrix.npz')
    
    if not os.path.exists(metadata_path):
        logger.error(f"Metadata not found: {metadata_path}")
        return None
    
    with open(metadata_path, 'r') as f:
        metadata = json.load(f)
    
    with open(gene_list_path, 'r') as f:
        gene_list = json.load(f)
    
    return {
        'metadata': metadata,
        'gene_list': gene_list,
        'expr_path': expr_path
    }

def compute_tsne_for_cells(expr_matrix, n_components=2, perplexity=30, 
                           learning_rate=200, max_iter=1000, random_state=42):
    """
    对表达矩阵计算t-SNE坐标
    """
    logger.info(f"Computing t-SNE for {expr_matrix.shape[0]} cells...")
    
    tsne = TSNE(
        n_components=n_components,
        perplexity=min(perplexity, expr_matrix.shape[0] - 1),
        learning_rate=learning_rate,
        max_iter=max_iter,
        random_state=random_state,
        metric='euclidean',
        init='pca',
        verbose=1
    )
    
    tsne_coords = tsne.fit_transform(expr_matrix)
    logger.info("t-SNE computation complete")
    
    return tsne_coords

def compute_tsne_for_subset(expr_matrix_sparse, cell_indices, gene_indices=None,
                            n_components=2, perplexity=30, max_iter=500):
    """
    对表达矩阵的子集计算t-SNE
    适用于筛选后的细胞子集
    """
    from scipy import sparse as sp
    
    n_cells = len(cell_indices)
    if n_cells < 10:
        logger.warning("Too few cells for t-SNE")
        return np.zeros((n_cells, n_components))
    
    # 提取子矩阵
    if gene_indices is not None:
        sub_matrix = expr_matrix_sparse[gene_indices, :][:, cell_indices].T
    else:
        sub_matrix = expr_matrix_sparse[:, cell_indices].T
    
    # 转为稠密矩阵
    if sp.issparse(sub_matrix):
        expr_dense = sub_matrix.toarray()
    else:
        expr_dense = sub_matrix
    
    # 调整perplexity
    actual_perplexity = min(perplexity, n_cells - 1)
    
    # 计算t-SNE
    tsne = TSNE(
        n_components=n_components,
        perplexity=actual_perplexity,
        max_iter=max_iter,
        random_state=42,
        init='pca',
        metric='euclidean'
    )
    
    tsne_coords = tsne.fit_transform(expr_dense)
    return tsne_coords

def get_tsne_coords_with_progress(expr_matrix, cell_indices=None, 
                                   callback=None, **tsne_params):
    """
    带进度的t-SNE计算
    """
    import time
    
    if cell_indices is not None:
        if hasattr(expr_matrix, 'toarray'):
            sub_data = expr_matrix[:, cell_indices].T.toarray()
        else:
            sub_data = expr_matrix[cell_indices, :]
    else:
        sub_data = expr_matrix
    
    n_cells = sub_data.shape[0]
    perplexity = min(tsne_params.get('perplexity', 30), n_cells - 1)
    max_iter = tsne_params.get('max_iter', 1000)
    
    start_time = time.time()
    
    tsne = TSNE(
        n_components=tsne_params.get('n_components', 2),
        perplexity=perplexity,
        max_iter=max_iter,
        learning_rate=tsne_params.get('learning_rate', 200),
        random_state=tsne_params.get('random_state', 42),
        init='pca',
        verbose=1
    )
    
    result = tsne.fit_transform(sub_data)
    
    elapsed = time.time() - start_time
    logger.info(f"t-SNE completed in {elapsed:.1f}s for {n_cells} cells")
    
    return result, elapsed

if __name__ == "__main__":
    print("Testing t-SNE module...")
    test_data = np.random.randn(100, 50)
    coords = compute_tsne_for_cells(test_data)
    print(f"Result shape: {coords.shape}")
    print("Test passed!")
