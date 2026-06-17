import anndata as ad

print("=== 检查Tabula Sapiens数据 ===")
adata = ad.read_h5ad("single_cell_data/tabula_sapiens_chunk0.h5ad")

print(f"细胞数: {adata.n_obs}")
print(f"基因数: {adata.n_vars}")

# 组织类型
tissues = adata.obs["organ_tissue"].unique().tolist()
print(f"\n组织类型 ({len(tissues)}个): {tissues[:15]}...")

# 细胞类型
print(f"细胞类型数: {adata.obs['cell_ontology_class'].nunique()}")

# 测试基因
test_genes = ["TP53", "BRCA1", "EGFR", "CD3D", "CD4", "CD8A", "MS4A1", "CD14", "LYZ", "GNLY", "NKG7", "KRT19", "EPCAM", "ALB", "ACTA2", "COL1A1"]

print(f"\n=== 测试基因覆盖 ===")
found = []
not_found = []
for gene in test_genes:
    if gene in adata.var_names:
        found.append(gene)
    else:
        matches = [g for g in adata.var_names if g.upper() == gene]
        if matches:
            found.append(f"{gene}({matches[0]})")
        else:
            not_found.append(gene)

print(f"找到: {found}")
print(f"未找到: {not_found}")

# 检查降维坐标
print(f"\n=== 降维坐标 ===")
for key in adata.obsm.keys():
    print(f"  {key}: {adata.obsm[key].shape}")
