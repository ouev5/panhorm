# PROJECT_MEMORY — 43 literature_search / GeneLit

## 项目定位
- 服务器：43.99.62.219
- 路径：`/www/wwwroot/literature_search`
- 服务：`literature-search.service`
- 端口：8005
- 现有功能：GeneLit 文献检索、DeepSeek 优化 PubMed 检索词、PubMed 文献检索、DeepSeek 提取基因/证据。

## AI 配置
- `app.py` 内置 DeepSeek 配置：
  - `DEEPSEEK_API_URL = https://api.deepseek.com/v1/chat/completions`
  - `DEEPSEEK_MODEL = deepseek-v4-pro`
- 注意：不要把 119 本机或别的 8005 当成这个服务；ouev5 要求改的是 43.99.62.219 的 8005。

## 2026-05-31 当前状态
- 新增 RAG PDF 报告能力：`POST /api/report/submit`
- 报告下载：`GET /api/report/download/{task_id}`
- 报告文件目录：`/www/wwwroot/literature_search/reports/`
- 结果 JSON 仍保存在：`/www/wwwroot/literature_search/results/`
- 已安装依赖：`reportlab`
- `uvicorn.run(... workers=2)` 已改为 `workers=1`，避免内存态 `tasks` 在多 worker 下轮询 404。

## 已验证示例
- 示例查询：`肢端肥大症的当前治疗策略`
- task_id：`a8f881f2`
- PDF：`/www/wwwroot/literature_search/reports/rag_report_a8f881f2.pdf`
- PDF 检查：22KB，PDF 1.4，2 pages
- SHA256：`5dc2981f12a3099b4e9674577296da1ab657bb6e85647b2245d9e75fc9cadd7c`
- DeepSeek 生成检索词：`acromegaly[MeSH] AND treatment[MeSH] AND strategy[Title/Abstract]`
- PubMed 候选文献数：12
- DeepSeek 报告标题：`肢端肥大症的当前治疗策略：多模式综合管理`

## 备份
- 改动前备份目录：`/www/wwwroot/literature_search/backups/`
- 主要备份：`app.py.before-rag-pdf-20260531-042123`
- 压缩备份：`literature_search_before_rag_pdf_20260531-042123.tar.gz`
