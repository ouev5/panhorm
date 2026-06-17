import scanpy as sc
import os

print("=== 查看Scanpy可用数据集 ===")

# 列出所有可用的数据集函数
dataset_funcs = [attr for attr in dir(sc.datasets) if not attr.startswith('_')]
print(f"可用数据集: {dataset_funcs}")

# 尝试下载更多数据
for name in dataset_funcs[:10]:
    print(f"\n=== 尝试 {name} ===")
    try:
        func = getattr(sc.datasets, name)
        adata = func()
        if hasattr(adata, 'n_obs'):
            print(f"细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
            adata.write_h5ad(f"single_cell_data/{name}.h5ad")
            size_mb = os.path.getsize(f"single_cell_data/{name}.h5ad") / 1024 / 1024
            print(f"保存成功: {size_mb:.1f} MB")
    except Exception as e:
        print(f"跳过: {type(e).__name__}")
