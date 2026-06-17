# 微调数据集使用说明

## 数据集概述

本数据集用于训练生物信息学分析平台（Gene Insights Hub + 动物激素数据库）的智能客服系统。

## 文件结构

```
fine_tuning_data/
├── project_info.json      # 项目基本信息
├── basic_qa.json          # 基础问答数据（10条）
├── advanced_qa.json       # 高级分析问答（8条）
├── troubleshooting.json   # 故障排除问答（8条）
├── hormone_qa.json        # 激素专业知识问答（8条）
├── workflow_qa.json       # 分析流程问答（8条）
└── rag_knowledge_base.json # RAG知识库文档（5篇）
```

## 数据统计

| 文件 | 问答对数量 | 用途 |
|------|-----------|------|
| basic_qa.json | 10 | 平台功能介绍 |
| advanced_qa.json | 8 | 结果解读指导 |
| troubleshooting.json | 8 | 问题解决 |
| hormone_qa.json | 8 | 激素专业知识 |
| workflow_qa.json | 8 | 操作流程指导 |
| **总计** | **42** | 微调训练 |

## 数据格式

### 问答数据格式（Alpaca格式）
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
      "content": "文档内容"
    }
  ]
}
```

## 使用建议

### 微调训练
1. 合并所有QA文件用于监督微调（SFT）
2. 建议使用Qwen2.5-7B或GLM-4作为基础模型
3. 训练轮数：3-5 epochs
4. 学习率：1e-5 到 5e-5

### RAG增强
1. 将rag_knowledge_base.json导入向量数据库
2. 推荐使用ChromaDB或FAISS
3. Embedding模型：BGE-large-zh或text-embedding-3-small
4. 检索top-k：3-5篇文档

### 系统架构建议
```
用户提问 → 向量检索(RAG) → 检索结果 + 问题 → 微调模型 → 回答
```

## 数据来源

所有数据基于以下项目代码和功能整理：
- 维恩图分析项目：/www/wwwroot/venn-tool/
- 动物激素数据库：/www/wwwroot/default/animal_hormone/

## 覆盖主题

1. **平台功能**（30%）
   - 维恩图分析
   - GO/KEGG富集
   - 单细胞分析
   - 桑基气泡图
   - STRING网络

2. **分析方法**（25%）
   - 结果解读
   - 参数选择
   - 最佳实践

3. **专业知识**（25%）
   - 激素分类
   - 信号通路
   - 疾病关联

4. **技术支持**（20%）
   - 故障排除
   - 使用技巧
   - 引用规范

## 扩展建议

可根据以下方向继续扩展数据集：
1. 添加更多具体疾病的案例分析
2. 增加不同分析场景的对话数据
3. 补充最新文献和方法学进展
4. 收集真实用户提问进行迭代

## 注意事项

1. 所有数据需确保科学准确性
2. 更新数据时保持格式一致
3. 定期验证回答的正确性
4. 注意数据库使用许可条款
