"""
生物信息学平台智能客服服务 - 完整版
RAG + DeepSeek API 混合架构
包含前端网页界面
"""

import os
import json
import requests
from typing import Optional, List, Dict
from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from sentence_transformers import SentenceTransformer
import numpy as np
from numpy import dot
from numpy.linalg import norm
import uvicorn

# ==================== 配置 ====================

# DeepSeek API配置
DEEPSEEK_API_KEY = "YOUR_VALUE_HERE"
DEEPSEEK_API_URL = "https://api.deepseek.com/v1/chat/completions"
DEEPSEEK_MODEL = "deepseek-v4-pro"

# 知识库路径
KNOWLEDGE_BASE_PATH = "/www/wwwroot/venn-tool/fine_tuning_data"

# 服务端口
SERVICE_PORT = 8002

# 系统提示词
SYSTEM_PROMPT = """You are an AI assistant for a bioinformatics analysis platform developed by Inner Mongolia University of Science and Technology.

The platform contains two main modules:
1. Gene Insights Hub（维恩图分析工具）- 支持 UpSet图、GO富集分析、KEGG分析、单细胞分析、桑基气泡图、STRING网络跳转
2. 动物激素调控数据库 - 激素数据查询、BLAST比对、转录组分析、DNA甲基化分析、GNN预测、系统发育分析、KEGG通路分析

Your responsibilities:
- 热情友好地回答用户关于平台功能的问题
- 解释生物信息学分析方法
- 提供激素相关知识（以当前 PanHorm 发布库与证据验证口径为准；旧库内容仅可作为候选线索，不能直接作为论文结论）
- 回答数据库规模/论文展示问题时，优先使用当前发布口径：126,199 records, 740 hormones, 1,978 genes, 3,716 diseases；激素→调控基因→基因对应疾病→DO标准化，激素-受体关系表单独维护。不要主动给出未复核旧数字。
- 指导用户操作步骤
- 解答常见问题

Response requirements:
- Professional, accurate, and friendly
- Answer in ENGLISH only
- Use clear structure with bullet points when appropriate
- Provide detailed steps for specific operations
"""

# 前端HTML
HTML_CONTENT = """
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>生物信息学平台智能客服</title>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh; display: flex; justify-content: center; align-items: center; padding: 20px;
        }
        .chat-container {
            width: 100%; max-width: 800px; background: #fff; border-radius: 20px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3); overflow: hidden;
            display: flex; flex-direction: column; height: 90vh; max-height: 700px;
        }
        .chat-header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white; padding: 20px; text-align: center;
        }
        .chat-header h1 { font-size: 1.5em; margin-bottom: 5px; }
        .chat-header p { font-size: 0.9em; opacity: 0.9; }
        .chat-messages {
            flex: 1; overflow-y: auto; padding: 20px; background: #f5f5f5;
        }
        .message { margin-bottom: 15px; display: flex; flex-direction: column; }
        .message.user { align-items: flex-end; }
        .message.assistant { align-items: flex-start; }
        .message-content {
            max-width: 80%; padding: 12px 18px; border-radius: 18px;
            line-height: 1.5; word-wrap: break-word; white-space: pre-wrap;
        }
        .message.user .message-content {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white; border-bottom-right-radius: 5px;
        }
        .message.assistant .message-content {
            background: white; color: #333; border-bottom-left-radius: 5px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        .chat-input {
            padding: 20px; background: white; border-top: 1px solid #eee;
            display: flex; gap: 10px;
        }
        .chat-input input {
            flex: 1; padding: 12px 20px; border: 2px solid #eee;
            border-radius: 25px; outline: none; font-size: 16px;
        }
        .chat-input input:focus { border-color: #667eea; }
        .chat-input button {
            padding: 12px 25px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white; border: none; border-radius: 25px;
            cursor: pointer; font-size: 16px;
        }
        .chat-input button:disabled { opacity: 0.6; cursor: not-allowed; }
        .typing {
            display: flex; gap: 5px; padding: 12px 18px; background: white;
            border-radius: 18px; box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        .typing span {
            width: 8px; height: 8px; background: #667eea; border-radius: 50%;
            animation: bounce 1.4s infinite ease-in-out;
        }
        .typing span:nth-child(1) { animation-delay: -0.32s; }
        .typing span:nth-child(2) { animation-delay: -0.16s; }
        @keyframes bounce {
            0%, 80%, 100% { transform: scale(0); }
            40% { transform: scale(1); }
        }
        .quick-questions {
            padding: 10px 20px; background: #f9f9f9;
            display: flex; flex-wrap: wrap; gap: 8px; border-top: 1px solid #eee;
        }
        .quick-question {
            padding: 8px 15px; background: white; border: 1px solid #ddd;
            border-radius: 20px; font-size: 13px; cursor: pointer;
        }
        .quick-question:hover {
            background: #667eea; color: white; border-color: #667eea;
        }
        @media (max-width: 600px) {
            .chat-container { height: 100vh; max-height: none; border-radius: 0; }
            .chat-header h1 { font-size: 1.2em; }
            .message-content { max-width: 90%; }
        }
    </style>
</head>
<body>
    <div class="chat-container">
        <div class="chat-header">
            <h1>🧬 生物信息学平台智能客服</h1>
            <p>Inner Mongolia University of Science and Technology</p>
        </div>
        <div class="quick-questions">
            <div class="quick-question" onclick="askQuestion('这个平台是做什么的？')">平台介绍</div>
            <div class="quick-question" onclick="askQuestion('如何进行GO富集分析？')">GO富集分析</div>
            <div class="quick-question" onclick="askQuestion('单细胞分析怎么用？')">单细胞分析</div>
            <div class="quick-question" onclick="askQuestion('BLAST比对怎么操作？')">BLAST比对</div>
        </div>
        <div class="chat-messages" id="chatMessages">
            <div class="message assistant">
                <div class="message-content">您好！我是 PanHorm / 生物信息学平台智能客服，可以帮您解答平台功能、分析方法和激素知识等问题。涉及数据库规模、激素-基因-疾病关系或论文展示时，我会优先按当前发布库与证据验证口径回答，旧库内容仅作为候选线索。请问有什么可以帮您的？😊</div>
            </div>
        </div>
        <div class="chat-input">
            <input type="text" id="userInput" placeholder="输入您的问题..." onkeypress="handleKeyPress(event)">
            <button id="sendBtn" onclick="sendMessage()">发送</button>
        </div>
    </div>
    <script>
        const chatMessages = document.getElementById('chatMessages');
        const userInput = document.getElementById('userInput');
        const sendBtn = document.getElementById('sendBtn');
        
        function addMessage(content, isUser) {
            const messageDiv = document.createElement('div');
            messageDiv.className = `message ${isUser ? 'user' : 'assistant'}`;
            messageDiv.innerHTML = `<div class="message-content">${content}</div>`;
            chatMessages.appendChild(messageDiv);
            chatMessages.scrollTop = chatMessages.scrollHeight;
        }
        
        function addTypingIndicator() {
            const typingDiv = document.createElement('div');
            typingDiv.className = 'message assistant';
            typingDiv.id = 'typingIndicator';
            typingDiv.innerHTML = '<div class="typing"><span></span><span></span><span></span></div>';
            chatMessages.appendChild(typingDiv);
            chatMessages.scrollTop = chatMessages.scrollHeight;
        }
        
        function removeTypingIndicator() {
            const typing = document.getElementById('typingIndicator');
            if (typing) typing.remove();
        }
        
        async function sendMessage() {
            const question = userInput.value.trim();
            if (!question) return;
            userInput.disabled = true;
            sendBtn.disabled = true;
            addMessage(question, true);
            userInput.value = '';
            addTypingIndicator();
            try {
                const response = await fetch('/api/chat', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ question: question })
                });
                const data = await response.json();
                removeTypingIndicator();
                addMessage(data.answer, false);
            } catch (error) {
                removeTypingIndicator();
                addMessage('抱歉，服务暂时出现问题，请稍后再试。', false);
            }
            userInput.disabled = false;
            sendBtn.disabled = false;
            userInput.focus();
        }
        
        function askQuestion(question) {
            userInput.value = question;
            sendMessage();
        }
        
        function handleKeyPress(event) {
            if (event.key === 'Enter') sendMessage();
        }
        userInput.focus();
    </script>
</body>
</html>
"""


# ==================== RAG检索模块 ====================

class RAGRetriever:
    def __init__(self, knowledge_base_path: str):
        print("正在加载embedding模型...")
        # 手动构建模型（解决本地加载Pooling参数缺失问题）
        from sentence_transformers import SentenceTransformer as ST, models
        _model_path = '/www/wwwroot/default/animal_hormone/local_models/all-MiniLM-L6-v2/'
        _transformer = models.Transformer(_model_path)
        _pooling = models.Pooling(
            _transformer.get_word_embedding_dimension(),
            pooling_mode_mean_tokens=True,
            pooling_mode_cls_token=False,
            pooling_mode_max_tokens=False
        )
        self.encoder = ST(modules=[_transformer, _pooling])
        self.knowledge_base_path = knowledge_base_path
        self.documents = []
        self.doc_embeddings = None
        self._load_knowledge_base()
        print(f"知识库加载完成，共 {len(self.documents)} 条知识")
    
    def _load_knowledge_base(self):
        kb_files = ["new_features_knowledge.json", "rag_knowledge_base.json", "rag_knowledge_base_extended.json"]
        for kb_file in kb_files:
            kb_path = os.path.join(self.knowledge_base_path, kb_file)
            if os.path.exists(kb_path):
                try:
                    with open(kb_path, 'r', encoding='utf-8') as f:
                        data = json.load(f)
                        if 'documents' in data:
                            for doc in data['documents']:
                                self.documents.append({
                                    'id': doc.get('id', ''),
                                    'title': doc.get('title', ''),
                                    'content': doc.get('content', '')
                                })
                    print(f"  已加载: {kb_file}")
                except: pass
        
        core_qa_files = ["basic_qa.json", "advanced_qa.json", "single_cell_qa.json",
            "hormone_qa.json", "operation_qa.json", "troubleshooting.json",
            "methodology.json", "missing_features_qa.json", "new_features_qa.json", "new_features_knowledge.json"]
        
        for qa_file in core_qa_files:
            qa_path = os.path.join(self.knowledge_base_path, qa_file)
            if os.path.exists(qa_path):
                try:
                    with open(qa_path, 'r', encoding='utf-8') as f:
                        qa_data = json.load(f)
                        if isinstance(qa_data, list):
                            for item in qa_data:
                                instruction = item.get('instruction', '')
                                output = item.get('output', '')
                                if instruction and output:
                                    self.documents.append({
                                        'id': instruction[:30],
                                        'title': instruction,
                                        'content': f"问：{instruction}\n答：{output}"
                                    })
                except: pass
        
        if self.documents:
            print(f"正在为 {len(self.documents)} 条知识生成向量...")
            texts = [doc['content'] for doc in self.documents]
            self.doc_embeddings = self.encoder.encode(texts, show_progress_bar=False)
            print("向量生成完成！")
    
    def search(self, query: str, top_k: int = 3) -> List[Dict]:
        if not self.documents or self.doc_embeddings is None: return []
        query_embedding = self.encoder.encode([query])
        similarities = []
        for i, doc_emb in enumerate(self.doc_embeddings):
            sim = dot(query_embedding[0], doc_emb) / (norm(query_embedding[0]) * norm(doc_emb))
            similarities.append((i, float(sim)))
        similarities.sort(key=lambda x: x[1], reverse=True)
        return [{**self.documents[i], 'score': round(sim, 4)} for i, sim in similarities[:top_k]]


def call_deepseek_api(messages: List[Dict], temperature: float = 0.7) -> str:
    headers = {"Content-Type": "application/json", "Authorization": f"Bearer {DEEPSEEK_API_KEY}"}
    data = {"model": DEEPSEEK_MODEL, "messages": messages, "temperature": temperature, "max_tokens": 2000}
    try:
        response = requests.post(DEEPSEEK_API_URL, headers=headers, json=data, timeout=30)
        response.raise_for_status()
        return response.json()["choices"][0]["message"]["content"]
    except: return "抱歉，服务暂时出现问题，请稍后再试。"


class CustomerService:
    def __init__(self):
        self.rag = RAGRetriever(KNOWLEDGE_BASE_PATH)
    
    def answer(self, question: str) -> Dict:
        relevant_docs = self.rag.search(question, top_k=3)
        context = ""
        if relevant_docs:
            context = "以下是与问题相关的参考资料：\n\n"
            for i, doc in enumerate(relevant_docs, 1):
                content = doc['content'][:500] + "..." if len(doc['content']) > 500 else doc['content']
                context += f"【参考{i}】{content}\n\n"
        
        messages = [{"role": "system", "content": SYSTEM_PROMPT}]
        if context:
            messages.append({"role": "system", "content": context + "\n请基于以上参考资料回答用户问题。"})
        messages.append({"role": "user", "content": question})
        
        return {
            "question": question,
            "answer": call_deepseek_api(messages),
            "sources": [{"title": d['title'][:50], "score": d['score']} for d in relevant_docs] if relevant_docs else []
        }


# ==================== FastAPI服务 ====================

app = FastAPI(title="生物信息学平台智能客服API", version="1.0")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

print("\n" + "="*50)
print("  初始化智能客服服务")
print("="*50 + "\n")
customer_service = CustomerService()
print("\n" + "="*50)
print("  服务就绪！")
print("="*50 + "\n")


class QuestionRequest(BaseModel):
    question: str

class AnswerResponse(BaseModel):
    question: str
    answer: str
    sources: List[Dict]


@app.get("/", response_class=HTMLResponse)
async def root():
    """返回聊天页面"""
    return HTML_CONTENT


@app.post("/api/chat", response_model=AnswerResponse)
async def chat(request: QuestionRequest):
    """对话接口"""
    if not request.question.strip():
        raise HTTPException(status_code=400, detail="问题不能为空")
    result = customer_service.answer(request.question)
    return AnswerResponse(**result)


@app.get("/api/health")
async def health():
    """健康检查"""
    return {
        "status": "ok",
        "service": "生物信息学平台智能客服",
        "knowledge_docs": len(customer_service.rag.documents),
        "api": "DeepSeek"
    }


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=SERVICE_PORT)
