# 微调数据集最终报告

## 项目概述

本数据集用于训练生物信息学分析平台智能客服系统，涵盖Gene Insights Hub维恩图分析工具和动物激素调控数据库两个项目。

## 数据集统计

### 总体数据量
- **问答对总数**: 501条
- **知识库文档**: 5篇
- **总文件大小**: 约280KB

### 问答数据分布

| 文件 | 数量 | 主题 |
|------|------|------|
| additional_qa.json | 45条 | 高级分析技巧、疾病/激素批量分析 |
| advanced_qa.json | 8条 | 结果解读指导 |
| basic_qa.json | 10条 | 平台功能介绍 |
| disease_analysis.json | 8条 | 疾病案例分析 |
| disease_gene_qa.json | 48条 | 疾病基因详细分析 |
| final_qa.json | 18条 | 操作步骤详解 |
| general_qa.json | 100条 | 通用知识问答 |
| hormone_qa_extended.json | 60条 | 激素专业知识扩展 |
| hormone_qa.json | 8条 | 激素基础知识 |
| last_batch_qa.json | 48条 | 最后补充问答 |
| methodology_deep_qa.json | 14条 | 方法学深入讲解 |
| methodology.json | 8条 | 方法学基础 |
| operation_qa.json | 25条 | 操作指南 |
| remaining_qa.json | 46条 | 剩余补充问答 |
| sankey_db_qa.json | 16条 | 桑基图/数据库操作 |
| single_cell_qa.json | 11条 | 单细胞分析 |
| troubleshooting_extended.json | 12条 | 故障排除扩展 |
| troubleshooting.json | 8条 | 常见问题 |
| workflow_qa.json | 8条 | 工作流程指导 |

### 知识库文档

| 文档ID | 标题 | 内容概述 |
|--------|------|----------|
| doc_001 | 平台功能总览 | 所有功能模块介绍 |
| doc_002 | 数据来源与可靠性 | 数据库引用信息 |
| doc_003 | 分析结果解读指南 | 结果解读方法 |
| doc_004 | 常见问题解决方案 | FAQ解答 |
| doc_005 | 激素数据库详细信息 | 数据库结构和内容 |

## 主题覆盖

### 1. 平台功能 (25%)
- 维恩图/UpSet图分析
- GO富集分析
- KEGG通路分析
- 桑基气泡图
- 单细胞聚类分析
- STRING蛋白互作网络
- 数据库搜索功能

### 2. 疾病研究 (20%)
- 20种常见疾病的基因分析流程
- 肿瘤研究方法
- 内分泌疾病
- 代谢疾病
- 免疫疾病
- 神经系统疾病

### 3. 激素知识 (15%)
- 20种重要激素详细解析
- 激素分类和功能
- 激素相关基因
- 激素信号通路

### 4. 方法学 (15%)
- 统计分析方法
- P值解读
- 多重检验校正
- 富集分析原理
- 降维算法

### 5. 操作指南 (15%)
- 详细操作步骤
- 使用技巧
- 快捷操作

### 6. 技术支持 (10%)
- 故障排除
- 常见问题解答
- 数据格式要求

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
top_k = 3
```

### 系统架构
```
用户提问
    ↓
意图识别
    ↓
向量检索(RAG) ← 知识库(5篇文档)
    ↓
检索结果 + 问题
    ↓
微调模型
    ↓
回答生成
    ↓
输出给用户
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

## 后续扩展建议

1. **对话数据**: 添加多轮对话模板
2. **工具调用**: 集成平台API调用能力
3. **实时更新**: 收集用户真实问答迭代
4. **多语言**: 添加英文问答数据
5. **专业知识**: 补充最新文献方法

---

生成时间: 2026-04-19
版本: 2.0
总计: 501条问答对 + 5篇知识库文档
