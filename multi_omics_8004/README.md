# autoBA_optimized - 生物信息学自动化分析平台

## 项目简介

autoBA_optimized 是一个基于 Python 3.13 的生物信息学自动化分析平台，旨在提供高效、智能的多组学数据分析能力。平台支持 RNA-seq、蛋白质组学、代谢组学和多组学数据的自动化分析，并集成了 LLM 技术，为用户提供智能的分析建议和结果解释。

### 核心功能

- **多组学数据分析**：支持 RNA-seq、蛋白质组学、代谢组学和多组学数据的分析
- **异步执行**：使用 asyncio 实现高效的任务并行处理
- **数据验证**：基于 Pydantic v2 进行严格的数据验证和类型检查
- **WebUI 界面**：使用 Gradio 提供友好的可视化操作界面
- **LLM 集成**：支持在线（OpenAI、Anthropic）和本地 LLM 模型
- **可配置的分析流程**：通过配置文件定义和管理分析流程

## 快速启动

### 环境要求

- Python 3.13
- pip 23.0+

### 安装步骤

1. **克隆项目**

   ```bash
   git clone <repository-url>
   cd autoBA_optimized
   ```

2. **创建虚拟环境**

   ```bash
   python -m venv venv
   # Windows
   venv\Scripts\activate
   # Linux/Mac
   source venv/bin/activate
   ```

3. **安装依赖**

   ```bash
   pip install -r requirements.txt
   ```

4. **配置 LLM**

   在 `config` 目录中创建 `llm_config.json` 文件，配置您的 LLM 提供商信息：

   ```json
   {
     "provider": "openai",
     "api_key": "YOUR_VALUE_HERE",
     "model": "gpt-4",
     "temperature": 0.7,
     "max_tokens": 1000
   }
   ```

5. **启动 WebUI**

   ```bash
   python src/webui/app.py
   ```

   然后在浏览器中访问 `http://localhost:7860` 开始使用平台。

## 目录结构

```
autoBA_optimized/
├── src/                 # 源代码目录
│   ├── __init__.py
│   ├── webui/           # WebUI 界面代码
│   │   └── __init__.py
├── config/              # 配置文件目录
│   └── pipelines/       # 分析流程配置
├── tests/               # 测试代码目录
│   └── fixtures/        # 测试数据
├── scripts/             # 辅助脚本
├── data/                # 数据存储目录
├── output/              # 分析结果输出目录
└── logs/                # 日志文件目录
```

## 分析流程

1. **数据上传**：通过 WebUI 上传您的组学数据文件
2. **选择分析类型**：选择 RNA-seq、蛋白质组学、代谢组学或多组学分析
3. **配置分析参数**：设置分析的具体参数和选项
4. **启动分析**：提交分析任务，系统会异步执行
5. **查看结果**：分析完成后，查看详细的分析结果和可视化图表
6. **获取智能建议**：基于 LLM 的分析建议和结果解释

## 技术栈

- **Python 3.13**：核心编程语言
- **asyncio**：异步执行框架
- **Pydantic v2**：数据验证和类型注解
- **Gradio**：WebUI 界面
- **OpenAI/Anthropic**：在线 LLM 支持
- **pandas/numpy**：数据处理
- **scikit-learn**：机器学习分析
- **matplotlib/seaborn**：数据可视化

## 贡献指南

欢迎对项目进行贡献！请按照以下步骤：

1. Fork 项目仓库
2. 创建您的特性分支 (`git checkout -b feature/amazing-feature`)
3. 提交您的更改 (`git commit -m 'Add some amazing feature'`)
4. 推送到分支 (`git push origin feature/amazing-feature`)
5. 打开一个 Pull Request

## 许可证

本项目采用 MIT 许可证。详见 LICENSE 文件。
