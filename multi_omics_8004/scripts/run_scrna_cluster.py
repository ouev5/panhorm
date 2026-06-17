#!/usr/bin/env python3

import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import os

# 解析命令行参数
parser = argparse.ArgumentParser(description='单细胞RNA-seq聚类分析')
parser.add_argument('--input', '-i', required=True, help='输入h5ad文件路径')
parser.add_argument('--sample-info', '-s', required=True, help='样本信息文件路径')
parser.add_argument('--output', '-o', required=True, help='输出h5ad文件路径')
parser.add_argument('--plot-dir', '-p', required=True, help='输出图片目录')
parser.add_argument('--n-top-genes', type=int, default=2000, help='用于PCA的top基因数')
parser.add_argument('--n-pcs', type=int, default=50, help='PCA主成分数')
parser.add_argument('--resolution', type=float, default=0.5, help='聚类分辨率')
parser.add_argument('--min-genes', type=int, default=200, help='每个细胞的最小基因数')
parser.add_argument('--min-cells', type=int, default=3, help='每个基因的最小细胞数')

args = parser.parse_args()

# 确保输出目录存在
os.makedirs(os.path.dirname(args.output), exist_ok=True)
os.makedirs(args.plot_dir, exist_ok=True)

# 读取数据
print("读取数据...")
adata = sc.read(args.input)

# 读取样本信息
sample_info = pd.read_csv(args.sample_info, index_col=0)

# 数据预处理
print("数据预处理...")
sc.pp.filter_cells(data, min_genes=args.min_genes)
sc.pp.filter_genes(data, min_cells=args.min_cells)

# 计算线粒体基因比例
mito_genes = data.var_names.str.startswith('MT-')
data.obs['percent_mito'] = sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=False)['percent_mt']

# 过滤线粒体基因比例高的细胞
data = data[data.obs['percent_mito'] < 5, :]

# 标准化和缩放
sc.pp.normalize_total(data, target_sum=1e4)
sc.pp.log1p(data)

# 选择高变基因
sc.pp.highly_variable_genes(data, n_top_genes=args.n_top_genes)
data = data[:, data.var.highly_variable]

# 批次校正（如果需要）
if 'batch' in data.obs:
    sc.pp.combat(data, key='batch')

# PCA降维
print("PCA降维...")
sc.tl.pca(data, n_comps=args.n_pcs)

# 计算邻接矩阵
sc.pp.neighbors(data, n_pcs=args.n_pcs)

# 聚类
print("聚类分析...")
sc.tl.leiden(data, resolution=args.resolution)

# UMAP可视化
print("UMAP可视化...")
sc.tl.umap(data)

# 保存UMAP图
plt.figure(figsize=(10, 8))
sc.pl.umap(data, color=['leiden'], show=False)
plt.savefig(os.path.join(args.plot_dir, 'umap_plot.png'), dpi=300, bbox_inches='tight')
plt.close()

# 计算每个聚类的标记基因
print("计算标记基因...")
sc.tl.rank_genes_groups(data, 'leiden', method='t-test')

# 保存标记基因
markers = sc.get.rank_genes_groups_df(data, group=None)
markers.to_csv(os.path.join(args.plot_dir, 'cluster_markers.csv'), index=False)

# 保存结果
print("保存结果...")
data.write(args.output)

print("分析完成！")
print(f"结果保存到: {args.output}")
print(f"UMAP图保存到: {os.path.join(args.plot_dir, 'umap_plot.png')}")
print(f"标记基因保存到: {os.path.join(args.plot_dir, 'cluster_markers.csv')}")
