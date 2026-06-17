import scanpy as sc
import os

print("=== Scanpy内置数据集 ===")
print("可用数据集:")
print("  - pbmc3k: 2700细胞")
print("  - pbmc68k_reduced: 700细胞(已处理)")
print("  - krasnow85: 肺部数据")
print("  - moignard15: 造血数据")
print("  - paul15: 造血数据")
print("  - burczynski19: 炎症数据")
print("")

# 尝试下载更多数据
datasets_to_try = [
    ("krasnow85", sc.datasets.krasnow85),
    ("moignard15", sc.datasets.moignard15),
    ("burczynski19", sc.datasets.burczynski19),
]

for name, func in datasets_to_try:
    print(f"=== 尝试下载 {name} ===")
    try:
        adata = func()
        print(f"细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
        adata.write_h5ad(f"single_cell_data/{name}.h5ad")
        size_mb = os.path.getsize(f"single_cell_data/{name}.h5ad") / 1024 / 1024
        print(f"保存成功: {size_mb:.1f} MB")
    except Exception as e:
        print(f"错误: {e}")
    print("")
