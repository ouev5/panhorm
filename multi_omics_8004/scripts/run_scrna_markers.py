#!/usr/bin/env python3

import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

# 解析命令行参数
parser = argparse.ArgumentParser(description='单细胞RNA-seq Marker基因鉴定')
parser.add_argument('--input', '-i', required=True, help='输入h5ad文件路径')
parser.add_argument('--cluster-column', '-c', required=True, help='聚类列名')
parser.add_argument('--output', '-o', required=True, help='输出marker基因文件路径')
parser.add_argument('--plot-dir', '-p', required=True, help='输出图片目录')
parser.add_argument('--method', default='wilcoxon', help='差异分析方法')
parser.add_argument('--n-genes', type=int, default=50, help='每个聚类的marker基因数')
parser.add_argument('--min-logfc', type=float, default=0.25, help='最小logFC')
parser.add_argument('--min-pvalue', type=float, default=0.05, help='最小p值')

args = parser.parse_args()

# 确保输出目录存在
os.makedirs(os.path.dirname(args.output), exist_ok=True)
os.makedirs(args.plot_dir, exist_ok=True)

# 读取数据
print("读取数据...")
adata = sc.read(args.input)

# 计算marker基因
print("计算marker基因...")
sc.tl.rank_genes_groups(adata, args.cluster_column, method=args.method)

# 获取marker基因
markers = sc.get.rank_genes_groups_df(adata, group=None)

# 过滤显著差异表达的基因
filtered_markers = markers[(markers['logfoldchanges'].abs() > args.min_logfc) & (markers['pvals_adj'] < args.min_pvalue)]

# 为每个聚类选择top N基因
top_markers = filtered_markers.groupby('group').head(args.n_genes)

# 保存marker基因
top_markers.to_csv(args.output, index=False)

# 生成热图
print("生成热图...")
top_genes_per_cluster = top_markers.groupby('group').head(10)['names'].tolist()
plt.figure(figsize=(15, 10))
sc.pl.heatmap(adata, var_names=top_genes_per_cluster, groupby=args.cluster_column, show=False)
plt.savefig(os.path.join(args.plot_dir, 'marker_heatmap.png'), dpi=300, bbox_inches='tight')
plt.close()

# 生成火山图（示例）
print("生成火山图...")
# 选择一个聚类作为示例
if len(adata.obs[args.cluster_column].unique()) > 0:
    example_cluster = adata.obs[args.cluster_column].unique()[0]
    cluster_markers = markers[markers['group'] == example_cluster]
    
    plt.figure(figsize=(10, 8))
    plt.scatter(cluster_markers['logfoldchanges'], -np.log10(cluster_markers['pvals_adj']), s=10)
    plt.xlabel('Log Fold Change')
    plt.ylabel('-Log10 p-value')
    plt.title(f'Marker Genes for Cluster {example_cluster}')
    plt.axhline(y=-np.log10(args.min_pvalue), color='r', linestyle='--')
    plt.axvline(x=args.min_logfc, color='r', linestyle='--')
    plt.axvline(x=-args.min_logfc, color='r', linestyle='--')
    plt.savefig(os.path.join(args.plot_dir, 'marker_volcano.png'), dpi=300, bbox_inches='tight')
    plt.close()

print("分析完成！")
print(f"Marker基因保存到: {args.output}")
print(f"热图保存到: {os.path.join(args.plot_dir, 'marker_heatmap.png')}")
print(f"火山图保存到: {os.path.join(args.plot_dir, 'marker_volcano.png')}")
