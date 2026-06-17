# 微调数据集最终报告 v4.0

## 项目概述

本数据集用于训练生物信息学分析平台智能客服系统，涵盖Gene Insights Hub维恩图分析工具和动物激素调控数据库两个项目。

**开发单位**：内蒙古科技大学生命科学与技术学院

## 数据集统计

### 总体数据量
- **问答对总数**: 806条
- **知识库文档**: 10篇
- **总文件大小**: 约450KB

### 问答数据分布

| 文件 | 数量 | 主题 |
|------|------|------|
| additional_qa.json | 45条 | 高级分析技巧、疾病/激素批量分析 |
| advanced_qa.json | 8条 | 结果解读指导 |
| basic_qa.json | 10条 | 平台功能介绍 |
| batch1-20_qa.json | 200条 | 新增补充问答（详见下表） |
| disease_analysis.json | 8条 | 疾病案例分析 |
| disease_gene_qa.json | 48条 | 疾病基因详细分析 |
| final_qa.json | 18条 | 操作步骤详解 |
| general_qa.json | 100条 | 通用知识问答 |
| hormone_qa_extended.json | 60条 | 激素专业知识扩展 |
| hormone_qa.json | 8条 | 激素基础知识 |
| last_batch_qa.json | 48条 | 最后补充问答 |
| methodology_deep_qa.json | 14条 | 方法学深入讲解 |
| methodology.json | 8条 | 方法学基础 |
| missing_features_qa.json | 25条 | BLAST/转录组/GNN/系统发育等 |
| operation_qa.json | 25条 | 操作指南 |
| remaining_qa.json | 46条 | 剩余补充问答 |
| sankey_db_qa.json | 16条 | 桑基图/数据库操作 |
| single_cell_qa.json | 11条 | 单细胞分析 |
| troubleshooting_extended.json | 12条 | 故障排除扩展 |
| troubleshooting.json | 8条 | 常见问题 |
| workflow_qa.json | 8条 | 工作流程指导 |

### 新增批次问答详情（batch1-20）

| 批次 | 数量 | 主题 |
|------|------|------|
| batch1 | 10条 | 实验室介绍、GO本体类型、基因集支持 |
| batch2 | 10条 | STRING跳转、单细胞数据来源、数据库字段 |
| batch3 | 10条 | 降维方法、基因表达查找、BLAST数据库 |
| batch4 | 10条 | 激素分类、胰岛素、皮质醇、甲状腺素等机制 |
| batch5 | 10条 | TSH、ACTH、促性腺激素、ADH、催产素等 |
| batch6 | 10条 | P值解读、多重检验校正、GO vs KEGG |
| batch7 | 10条 | 单细胞质控、聚类、细胞类型识别 |
| batch8 | 10条 | BLAST程序选择、转录组重复、差异基因 |
| batch9 | 10条 | 可视化方式、图表下载、统计方法选择 |
| batch10 | 10条 | 物种支持、基因名称处理、API接口 |
| batch11 | 10条 | 基因表达、转录因子、启动子、增强子 |
| batch12 | 10条 | UniProt、PubMed、PMID、Ensembl、HGNC |
| batch13 | 10条 | 激素受体分类、GPCR、核受体、激素抵抗 |
| batch14 | 10条 | 疾病分析：糖尿病、肿瘤、心血管、PCOS等 |
| batch15 | 10条 | Enrichment Score、标准化、批次效应、GSEA |
| batch16 | 10条 | 翻译后修饰、多态性、外显子、miRNA、lncRNA |
| batch17 | 10条 | 平台协同使用、课题设计、引用方式 |
| batch18 | 10条 | 蛋白质功能预测、信号通路、多组学整合 |
| batch19 | 10条 | 学校介绍、数据来源、商业使用、论文描述 |
| batch20 | 10条 | NES、leading edge、时空特异性、基因印记 |

### 知识库文档

| 文档ID | 标题 | 内容概述 |
|--------|------|----------|
| doc_001 | 平台功能总览 | 所有功能模块介绍 |
| doc_002 | 数据来源与可靠性 | 数据库引用信息 |
| doc_003 | 分析结果解读指南 | 结果解读方法 |
| doc_004 | 常见问题解决方案 | FAQ解答 |
| doc_005 | 激素数据库详细信息 | 数据库结构和内容 |
| doc_006 | BLAST比对功能详细说明 | BLAST程序类型、参数、结果解读 |
| doc_007 | 转录组分析功能详解 | 分析流程、可视化图表、结果下载 |
| doc_008 | GNN基因网络预测功能说明 | 预测原理、置信度解读、应用场景 |
| doc_009 | 系统发育与跨物种分析功能 | 进化树构建、保守性分析 |
| doc_010 | AI功能详解：RAG检索与论文问答 | AI系统架构、使用方法、对比说明 |

## 主题覆盖

### 1. 平台与实验室介绍 (5%)
- 内蒙古科技大学生命科学与技术学院介绍
- 平台开发目的和使用说明
- 引用和联系方式

### 2. 激素专业知识 (20%)
- 激素分类（肽类、类固醇、胺类）
- 各种激素详细机制
- 激素受体类型和信号通路
- 激素相关疾病

### 3. 富集分析 (15%)
- GO富集分析原理和解读
- KEGG通路分析
- P值和多重检验
- 结果筛选和验证

### 4. 单细胞分析 (10%)
- 数据质控
- 降维方法选择
- 细胞类型识别
- 差异表达分析

### 5. 其他分析功能 (20%)
- BLAST比对
- 转录组分析
- DNA甲基化
- GNN预测
- 系统发育分析
- AI功能

### 6. 数据库知识 (10%)
- UniProt、PubMed、KEGG等
- 基因ID类型
- 数据引用规范

### 7. 方法学知识 (10%)
- 统计分析原理
- 实验设计
- 数据处理方法

### 8. 应用指南 (10%)
- 疾病研究应用
- 论文写作指南
- 问题解决方案

## 数据格式

### 问答数据 (Alpaca格式)
```json
{
  "instruction": "问题文本",
  "input": "",
  "output": "回答文本"
}
```

### RAG知识库格式
```json
{
  "documents": [
    {
      "id": "doc_xxx",
      "title": "文档标题",
      "content": "详细内容"
    }
  ]
}
```

## 使用建议

### 微调训练配置
```python
# 推荐配置
base_model = "Qwen/Qwen2.5-7B-Instruct"
# 或
base_model = "THUDM/glm-4-9b-chat"

# 训练参数
epochs = 3-5
learning_rate = 2e-5
batch_size = 8
warmup_ratio = 0.1
```

### RAG系统配置
```python
# 向量数据库
vector_db = "ChromaDB"  # 或 FAISS

# Embedding模型
embedding_model = "BAAI/bge-large-zh-v1.5"

# 检索参数
chunk_size = 500
chunk_overlap = 50
top_k = 5
```

## 文件位置

- **服务器路径**: `/www/wwwroot/venn-tool/fine_tuning_data/`
- **服务器IP**: 43.99.62.219

## 数据来源

所有数据基于以下项目整理：
- 维恩图分析项目: `/www/wwwroot/venn-tool/`
- 动物激素数据库: `/www/wwwroot/default/animal_hormone/`
- 数据库: MySQL `animal_hormone` (4,096条激素记录)
- 单细胞数据: Tabula Sapiens (94,836细胞)

---

生成时间: 2026-04-19
版本: 4.0
总计: 806条问答对 + 10篇知识库文档
开发单位: 内蒙古科技大学生命科学与技术学院
