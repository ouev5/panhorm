#!/bin/bash
# 快速安装生物信息学工具
# 用法：bash install_bio.sh

echo "=== 快速安装生物信息学工具 ==="

# 设置清华镜像加速
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
conda config --set show_channel_urls yes

# 一次性安装所有包
echo "正在安装所有工具，请稍候..."
conda install -y -c bioconda \
    samtools \
    fastqc \
    hisat2 \
    subread \
    multiqc

conda install -y -c conda-forge \
    r-base \
    r-ggplot2 \
    r-dplyr \
    r-tidyr \
    bioconductor-deseq2

conda install -y \
    pandas \
    numpy \
    matplotlib \
    seaborn \
    scikit-learn

echo "=== 安装完成 ==="