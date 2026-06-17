from django.db.models import Q
from hormone_app.models import HormoneRelatedGene, HormoneReceptorFull, HormoneReceptorInfo
import numpy as np
from sentence_transformers import SentenceTransformer
import chromadb
from chromadb.config import Settings
import threading
import re
from itertools import chain
from functools import lru_cache
from datetime import datetime
import hashlib
import json
import requests
import logging
import sys
import os
import traceback

# ==================== 配置日志 ====================
# 创建日志目录
log_dir = '/www/wwwroot/default/animal_hormone/logs'
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

# 配置日志 - 修正 FileHandler 的用法
logger = logging.getLogger('rag_deepseek')
logger.setLevel(logging.INFO)

# 控制台处理器
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
console_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(console_format)
logger.addHandler(console_handler)

# 文件处理器 - 普通日志
file_handler = logging.FileHandler(os.path.join(log_dir, 'rag_deepseek.log'), encoding='utf-8')
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(console_format)
logger.addHandler(file_handler)

# 错误日志处理器
error_file_handler = logging.FileHandler(os.path.join(log_dir, 'rag_deepseek_error.log'), encoding='utf-8')
error_file_handler.setLevel(logging.ERROR)
error_file_handler.setFormatter(console_format)
logger.addHandler(error_file_handler)

# ==================== 配置常量 ====================
FIELD_MAPPING = {
    "hormone_name": "Hormone Name",
    "related_genes": "Related Genes",
    "pmid": "PMID",
    "gene_sequence": "Gene Sequence",
    "related_diseases": "Related Diseases",
    "hormone_uniprot_id": "Hormone UniProt ID",
    "hormone_species_name": "Hormone Species Name",
    "receptor_name": "Receptor Name",
    "receptor_uniprot_id": "Receptor UniProt ID",
    "pubchem_id": "PubChem ID",
    "receptor_coding_genes": "Receptor Coding Genes",
    "protein_sequence": "Protein Sequence",
    "dna_sequence": "DNA Sequence",
    "mrna_sequence": "mRNA Sequence",
    "amino_acid_sequence": "Amino Acid Sequence",
    "receptor_coding_genes_sequence": "Receptor Coding Genes Sequence",
    "hormone_coding_genes_sequence": "Hormone Coding Genes Sequence"
}

SEQUENCE_FIELDS = [
    "gene_sequence", "protein_sequence", "dna_sequence", 
    "mrna_sequence", "amino_acid_sequence", "sequence",
    "receptor_coding_genes_sequence", "hormone_coding_genes_sequence"
]

# Serper 搜索 API 配置
SERPER_API_KEY = "YOUR_VALUE_HERE"
SERPER_API_URL = "google.serper.dev"
SERPER_API_PATH = "/search"

CHUNK_SIZE = 512
CHUNK_OVERLAP = 50
CACHE_TTL = 300

# DeepSeek配置
DEEPSEEK_API_KEY = "YOUR_VALUE_HERE"
DEEPSEEK_API_URL = "https://api.deepseek.com/v1/chat/completions"
DEEPSEEK_MODEL = "deepseek-v4-pro"

class DeepSeekRAGClient:
    """DeepSeek RAG客户端 - 整合检索和生成"""
    
    def __init__(self, api_key=DEEPSEEK_API_KEY, model=DEEPSEEK_MODEL):
        self.api_key = api_key
        self.model = model
        self.session = requests.Session()
        self.headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {self.api_key}"
        }
        logger.info(f"DeepSeekRAGClient initialized with model: {model}")
    
    def generate_with_context(self, question, context, context_type="semantic"):
        """使用上下文生成回答"""
        logger.info(f"Generating answer with DeepSeek - context_type: {context_type}")
        logger.info(f"Question: {question}")
        logger.info(f"Context length: {len(context)} characters")
        logger.debug(f"Context preview: {context[:500]}...")
        
        # 系统提示
        system_prompts = {
            "structured": "You are a hormone database assistant. Provide concise and accurate information based on the database records. If no relevant information is found, state that clearly.",
            "semantic": "You are a hormone research assistant. Answer questions based on the provided context from scientific literature and database records. Use only the provided context, do not add external knowledge.",
            "extracted": "You are a hormone information assistant. Use the extracted information to answer the question. If the information is incomplete, explain what is missing.",
            "web_search": "You are a helpful assistant. Use the provided web search results to answer the question. If the results are insufficient or not directly relevant, say so. Do not make up information.",
            "hybrid": "You are a hormone research assistant. The following context contains both database records and web search results. Use all provided information to answer the question. If possible, distinguish between database-sourced and web-sourced information in your answer. If the information is insufficient, explain what is missing."
        }
        
        system_prompt = system_prompts.get(context_type, "You are a professional hormone research assistant. Answer questions based on the provided context.")
        logger.debug(f"System prompt: {system_prompt}")
        
        # 构建消息
        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": f"""Based on the following context, please answer the question:

Context:
{context}

Question: {question}

Guidelines:
- Use ONLY the information provided in the context
- If the information is insufficient, state what's missing
- Be concise but informative
- Answer in English
"""}
        ]
        
        # API参数
        api_params = {
            "model": self.model,
            "messages": messages,
            "max_tokens": 512,
            "temperature": 0.3,
            "top_p": 0.95,
            "stream": False
        }
        
        logger.info(f"Sending request to DeepSeek API, model: {self.model}")
        
        try:
            response = self.session.post(
                url=DEEPSEEK_API_URL,
                json=api_params,
                headers=self.headers,
                timeout=30,
                verify=True
            )
            
            logger.info(f"DeepSeek API response status: {response.status_code}")
            
            if response.status_code == 200:
                api_result = response.json()
                if "choices" in api_result and len(api_result["choices"]) > 0:
                    reply = api_result["choices"][0]["message"]["content"].strip()
                    usage = api_result.get("usage", {})
                    logger.info(f"DeepSeek generation successful, response length: {len(reply)} characters")
                    logger.info(f"Token usage: {usage}")
                    logger.debug(f"Response preview: {reply[:200]}...")
                    
                    return {
                        "success": True,
                        "response": reply,
                        "model": self.model,
                        "usage": usage
                    }
                else:
                    error_msg = "API returned no choices"
                    logger.error(f"DeepSeek API error: {error_msg}")
                    logger.error(f"API response: {api_result}")
                    return {"success": False, "error": error_msg}
            else:
                error_msg = self._handle_error(response)
                logger.error(f"DeepSeek API error: {error_msg}")
                return {"success": False, "error": error_msg}
            
        except requests.exceptions.Timeout:
            logger.error("DeepSeek API timeout")
            return {"success": False, "error": "Request timeout"}
        except requests.exceptions.ConnectionError:
            logger.error("DeepSeek API connection error")
            return {"success": False, "error": "Connection error"}
        except Exception as e:
            logger.error(f"DeepSeek API unexpected error: {str(e)}")
            logger.error(traceback.format_exc())
            return {"success": False, "error": str(e)}
    
    def _handle_error(self, response):
        """处理API错误"""
        status_code = response.status_code
        try:
            error_detail = response.json()
            error_message = error_detail.get("error", {}).get("message", response.text)
        except:
            error_message = response.text[:200]
        
        if status_code == 401:
            return "API Key is invalid or expired"
        elif status_code == 429:
            return "Rate limit exceeded or quota exhausted"
        elif status_code == 404:
            return "API endpoint error"
        else:
            return f"API error (status {status_code}): {error_message}"

# ==================== ChromaDB向量存储 ====================
class HormoneChromaStore:
    """基于ChromaDB的向量存储"""
    
    def __init__(self, collection_name="hormone_knowledge"):
        logger.info(f"Initializing HormoneChromaStore with collection: {collection_name}")
        
        self.client = chromadb.Client(Settings(
            persist_directory="./chromadb_data",
            anonymized_telemetry=False
        ))
        
        try:
            self.collection = self.client.get_collection(collection_name)
            count = self.collection.count()
            logger.info(f"Loaded existing ChromaDB collection: {collection_name}, count: {count}")
        except Exception as e:
            logger.info(f"Collection {collection_name} not found, creating new one")
            self.collection = self.client.create_collection(
                name=collection_name,
                metadata={"hnsw:space": "cosine"}
            )
            logger.info(f"Created new ChromaDB collection: {collection_name}")
        
        logger.info("Loading SentenceTransformer model: BAAI/bge-small-en-v1.5")
        self.model = SentenceTransformer('BAAI/bge-small-en-v1.5')
        logger.info("Model loaded successfully")
        
        self.query_cache = {}
        self.cache_timestamps = {}
        logger.info("HormoneChromaStore initialization complete")
    
    def _split_text_into_chunks(self, text, metadata):
        """将长文本分割成块"""
        if not text or len(text) < CHUNK_SIZE:
            return [(text, metadata)]
        
        words = text.split()
        chunks = []
        
        for i in range(0, len(words), CHUNK_SIZE - CHUNK_OVERLAP):
            chunk_words = words[i:i + CHUNK_SIZE]
            chunk_text = ' '.join(chunk_words)
            
            chunk_metadata = metadata.copy()
            chunk_metadata['chunk_index'] = i // (CHUNK_SIZE - CHUNK_OVERLAP)
            chunk_metadata['total_chunks'] = (len(words) + CHUNK_SIZE - 1) // CHUNK_SIZE
            
            chunks.append((chunk_text, chunk_metadata))
        
        logger.debug(f"Split text into {len(chunks)} chunks")
        return chunks
    
    def rebuild_index(self):
        """从数据库重建向量索引"""
        logger.info("=" * 50)
        logger.info("Starting to rebuild ChromaDB vector index...")
        
        try:
            documents = []
            metadatas = []
            ids = []
            
            # HormoneRelatedGene table
            logger.info("Processing HormoneRelatedGene table...")
            gene_count = 0
            for item in HormoneRelatedGene.objects.all():
                if item.hormone_name:
                    if item.related_genes and item.related_genes != "Not available":
                        text = f"Hormone {item.hormone_name} is associated with genes: {item.related_genes}"
                        chunks = self._split_text_into_chunks(text, {
                            'source': 'HormoneRelatedGene',
                            'db_id': item.id,
                            'hormone': item.hormone_name,
                            'field': 'related_genes'
                        })
                        
                        for chunk_text, chunk_meta in chunks:
                            doc_id = f"hrg_gene_{item.id}_{chunk_meta.get('chunk_index', 0)}"
                            documents.append(chunk_text)
                            metadatas.append(chunk_meta)
                            ids.append(doc_id)
                            gene_count += 1
                    
                    if item.related_diseases and item.related_diseases != "Not available":
                        text = f"Hormone {item.hormone_name} is related to diseases: {item.related_diseases}"
                        chunks = self._split_text_into_chunks(text, {
                            'source': 'HormoneRelatedGene',
                            'db_id': item.id,
                            'hormone': item.hormone_name,
                            'field': 'related_diseases'
                        })
                        
                        for chunk_text, chunk_meta in chunks:
                            doc_id = f"hrg_disease_{item.id}_{chunk_meta.get('chunk_index', 0)}"
                            documents.append(chunk_text)
                            metadatas.append(chunk_meta)
                            ids.append(doc_id)
                            gene_count += 1
            
            logger.info(f"HormoneRelatedGene: {gene_count} chunks generated")
            
            # HormoneReceptorFull table
            logger.info("Processing HormoneReceptorFull table...")
            receptor_full_count = 0
            for item in HormoneReceptorFull.objects.all():
                if item.hormone_name and item.receptor_name:
                    text = f"Hormone {item.hormone_name} binds to receptor {item.receptor_name}"
                    chunks = self._split_text_into_chunks(text, {
                        'source': 'HormoneReceptorFull',
                        'db_id': item.id,
                        'hormone': item.hormone_name,
                        'receptor': item.receptor_name,
                        'field': 'receptor_binding'
                    })
                    
                    for chunk_text, chunk_meta in chunks:
                        doc_id = f"hrf_{item.id}_{chunk_meta.get('chunk_index', 0)}"
                        documents.append(chunk_text)
                        metadatas.append(chunk_meta)
                        ids.append(doc_id)
                        receptor_full_count += 1
            
            logger.info(f"HormoneReceptorFull: {receptor_full_count} chunks generated")
            
            # HormoneReceptorInfo table
            logger.info("Processing HormoneReceptorInfo table...")
            receptor_info_count = 0
            for item in HormoneReceptorInfo.objects.all():
                if item.hormone_name and item.receptor_name:
                    text = f"Hormone {item.hormone_name} interacts with receptor {item.receptor_name}"
                    if item.pubchem_id and item.pubchem_id != "Not available":
                        text += f" (PubChem ID: {item.pubchem_id})"
                    
                    chunks = self._split_text_into_chunks(text, {
                        'source': 'HormoneReceptorInfo',
                        'db_id': item.id,
                        'hormone': item.hormone_name,
                        'receptor': item.receptor_name,
                        'field': 'receptor_interaction'
                    })
                    
                    for chunk_text, chunk_meta in chunks:
                        doc_id = f"hri_{item.id}_{chunk_meta.get('chunk_index', 0)}"
                        documents.append(chunk_text)
                        metadatas.append(chunk_meta)
                        ids.append(doc_id)
                        receptor_info_count += 1
            
            logger.info(f"HormoneReceptorInfo: {receptor_info_count} chunks generated")
            
            if not documents:
                logger.warning("No data for vector index")
                return False
            
            total_chunks = len(documents)
            logger.info(f"Total chunks to process: {total_chunks}")
            logger.info(f"Generating embeddings...")
            
            batch_size = 64
            for i in range(0, len(documents), batch_size):
                batch_end = min(i + batch_size, len(documents))
                batch_docs = documents[i:batch_end]
                batch_ids = ids[i:batch_end]
                batch_metas = metadatas[i:batch_end]
                
                logger.info(f"Processing batch {i//batch_size + 1}/{(total_chunks + batch_size - 1)//batch_size}")
                
                embeddings = self.model.encode(batch_docs).tolist()
                
                self.collection.add(
                    embeddings=embeddings,
                    documents=batch_docs,
                    metadatas=batch_metas,
                    ids=batch_ids
                )
                
                logger.info(f"Processed {batch_end}/{total_chunks} chunks")
            
            final_count = self.collection.count()
            logger.info(f"Vector index built successfully with {final_count} chunks")
            logger.info("=" * 50)
            return True
            
        except Exception as e:
            logger.error(f"Error rebuilding index: {str(e)}")
            logger.error(traceback.format_exc())
            return False
    
    def semantic_search(self, query: str, top_k: int = 5):
        """语义搜索"""
        logger.info(f"Executing semantic search - query: '{query}', top_k: {top_k}")
        
        try:
            # 检查缓存
            cache_key = hashlib.md5(f"{query}_{top_k}".encode()).hexdigest()
            if cache_key in self.query_cache:
                cache_time = self.cache_timestamps.get(cache_key)
                if cache_time and (datetime.now() - cache_time).seconds < CACHE_TTL:
                    logger.info(f"Using cached results for: '{query}', cache age: {(datetime.now() - cache_time).seconds}s")
                    cached_results = self.query_cache[cache_key]
                    logger.info(f"Cache returned {len(cached_results)} results")
                    return cached_results
            
            logger.info("Generating query embedding...")
            query_embedding = self.model.encode([query]).tolist()
            logger.info(f"Query embedding generated, dimension: {len(query_embedding[0])}")
            
            collection_count = self.collection.count()
            logger.info(f"Querying ChromaDB collection (size: {collection_count})...")
            
            results = self.collection.query(
                query_embeddings=query_embedding,
                n_results=min(top_k * 2, collection_count),
                include=["documents", "metadatas", "distances"]
            )
            
            logger.info(f"ChromaDB returned {len(results['documents'][0]) if results['documents'] else 0} raw results")
            
            processed_results = []
            seen_records = set()
            
            if results['documents'] and results['documents'][0]:
                for i, (doc, metadata, distance) in enumerate(zip(
                    results['documents'][0],
                    results['metadatas'][0],
                    results['distances'][0]
                )):
                    similarity = float(1 / (1 + distance))
                    record_key = f"{metadata['source']}_{metadata['db_id']}"
                    
                    logger.debug(f"Result {i+1}: distance={distance:.4f}, similarity={similarity:.4f}, source={metadata['source']}, db_id={metadata['db_id']}")
                    logger.debug(f"Text preview: {doc[:100]}...")
                    
                    if record_key not in seen_records:
                        seen_records.add(record_key)
                        processed_results.append({
                            'text': doc,
                            'similarity': similarity,
                            'metadata': metadata
                        })
                        logger.debug(f"Added result {len(processed_results)}")
                        
                        if len(processed_results) >= top_k:
                            logger.info(f"Reached target of {top_k} results, stopping")
                            break
            
            logger.info(f"Semantic search completed, returning {len(processed_results)} results")
            
            # 更新缓存
            self.query_cache[cache_key] = processed_results
            self.cache_timestamps[cache_key] = datetime.now()
            logger.info(f"Results cached with key: {cache_key}")
            
            return processed_results
            
        except Exception as e:
            logger.error(f"Semantic search error: {str(e)}")
            logger.error(traceback.format_exc())
            return []

def serper_search(query: str, num_results: int = 5):
    """使用 Serper API 进行网络搜索，返回格式化的文本结果"""
    logger.info(f"Performing Serper web search for: '{query}'")
    
    import http.client
    import json
    
    payload = json.dumps({
        "q": query,
        "num": num_results  # 请求返回的结果数量
    })
    
    headers = {
        'X-API-KEY': SERPER_API_KEY,
        'Content-Type': 'application/json'
    }
    
    try:
        conn = http.client.HTTPSConnection(SERPER_API_URL, timeout=10)
        conn.request("POST", SERPER_API_PATH, payload, headers)
        res = conn.getresponse()
        data = res.read()
        
        if res.status != 200:
            logger.error(f"Serper API error: status {res.status}, response: {data.decode('utf-8')}")
            return None
        
        result_json = json.loads(data.decode("utf-8"))
        logger.info(f"Serper API response parsed successfully")
        
        # 提取 organic 结果（自然搜索结果）
        organic = result_json.get("organic", [])
        logger.info(f"Serper returned {len(organic)} organic results")
        
        if not organic:
            logger.info("Serper returned no organic results")
            return None  # 返回None而不是字符串，方便判断
        
        # 构建格式化的文本
        lines = []
        for i, item in enumerate(organic[:num_results], 1):
            title = item.get("title", "No title")
            link = item.get("link", "")
            snippet = item.get("snippet", "")
            lines.append(f"{i}. {title}\n   URL: {link}\n   Summary: {snippet}\n")
        
        formatted_results = "\n".join(lines)
        logger.info(f"Formatted {len(lines)} search results")
        return formatted_results
        
    except Exception as e:
        logger.error(f"Serper search failed: {str(e)}")
        logger.error(traceback.format_exc())
        return None

# ==================== 关键词检索 ====================
def keyword_retrieval(query: str):
    """关键词检索"""
    logger.info(f"Executing keyword retrieval for query: '{query}'")
    
    if not query or query.strip() == "":
        logger.warning("Empty query provided")
        return None

    def format_item(item):
        formatted = {}
        for k, v in item.items():
            field_name = FIELD_MAPPING.get(k, k)
            if any(seq_field in k.lower() for seq_field in SEQUENCE_FIELDS) or "sequence" in k.lower():
                if v and isinstance(v, str) and len(v) > 0:
                    formatted[field_name] = "Sequence data is hidden for security"
                else:
                    formatted[field_name] = "Not available"
            else:
                formatted[field_name] = v if v is not None else "Not available"
        return formatted

    # HormoneRelatedGene table
    hgd_query = HormoneRelatedGene.objects.filter(
        Q(hormone_name__iexact=query) |
        Q(hormone_name__icontains=query) |
        Q(related_genes__icontains=query) |
        Q(related_diseases__icontains=query)
    )
    hgd_data = [format_item(item) for item in hgd_query.values()] if hgd_query.exists() else []
    logger.info(f"HormoneRelatedGene: {len(hgd_data)} results")

    # HormoneReceptorFull table
    phr_query = HormoneReceptorFull.objects.filter(
        Q(hormone_name__iexact=query) |
        Q(hormone_name__icontains=query) |
        Q(receptor_name__iexact=query) |
        Q(receptor_name__icontains=query)
    )
    phr_data = [format_item(item) for item in phr_query.values()] if phr_query.exists() else []
    logger.info(f"HormoneReceptorFull: {len(phr_data)} results")

    # HormoneReceptorInfo table
    nphr_query = HormoneReceptorInfo.objects.filter(
        Q(hormone_name__iexact=query) |
        Q(hormone_name__icontains=query) |
        Q(receptor_name__iexact=query) |
        Q(receptor_name__icontains=query)
    )
    nphr_data = [format_item(item) for item in nphr_query.values()] if nphr_query.exists() else []
    logger.info(f"HormoneReceptorInfo: {len(nphr_data)} results")

    result = {
        "hormone_related_gene": hgd_data,
        "hormone_receptor_full": phr_data,
        "hormone_receptor_info": nphr_data
    }
    
    total_results = len(hgd_data) + len(phr_data) + len(nphr_data)
    logger.info(f"Keyword retrieval total results: {total_results}")
    
    return result

# ==================== 初始化 ====================
logger.info("=" * 50)
logger.info("Initializing RAG DeepSeek module...")
vector_store = HormoneChromaStore()
deepseek_rag_client = DeepSeekRAGClient()
logger.info("Module initialization complete")
logger.info("=" * 50)

def initialize_vector_store_async():
    """异步初始化向量存储"""
    logger.info("Starting background thread for vector index loading...")
    
    def load_in_background():
        try:
            logger.info("Background thread: checking vector store...")
            count = vector_store.collection.count()
            logger.info(f"Current vector store count: {count}")
            
            if count == 0:
                logger.info("Vector store empty, starting rebuild...")
                vector_store.rebuild_index()
            else:
                logger.info("Vector store already has data, skipping rebuild")
                
        except Exception as e:
            logger.error(f"Background loading failed: {str(e)}")
            logger.error(traceback.format_exc())
    
    thread = threading.Thread(target=load_in_background, daemon=True)
    thread.start()
    logger.info("Background thread started")

# ==================== 查询分析 ====================
def is_simple_query(query: str) -> bool:
    """判断是否为简单查询（直接查找激素/基因）"""
    query = query.strip()
    if len(query) == 0:
        return False
    
    words = query.split()
    if len(words) > 3:
        logger.debug(f"Query has {len(words)} words > 3, considered complex")
        return False
    
    question_words = {'which', 'what', 'how', 'why', 'where', 'when', 
                     'is', 'are', 'does', 'do', 'did', 'can', 'could'}
    if any(w.lower() in question_words for w in words):
        logger.debug(f"Query contains question word, considered complex")
        return False
    
    if '?' in query:
        logger.debug(f"Query contains question mark, considered complex")
        return False
    
    logger.debug(f"Query considered simple")
    return True

# ==================== 上下文构建 ====================
def build_context_from_keyword(keyword_results: dict):
    """从关键词结果构建上下文"""
    logger.info("Building context from keyword results...")
    
    lines = []
    for table_name, records in keyword_results.items():
        for rec in records:
            valid_fields = [
                f"{k}: {v}" 
                for k, v in rec.items() 
                if v != "Not available" and "Sequence data is hidden" not in str(v)
            ]
            if valid_fields:
                lines.append(" | ".join(valid_fields))
    
    context = "\n".join(lines) if lines else None
    logger.info(f"Keyword context built, length: {len(context) if context else 0} characters")
    if context:
        logger.debug(f"Keyword context preview: {context[:200]}...")
    
    return context

def build_context_from_semantic(semantic_results: list, keyword_results: dict = None):
    """从语义结果构建上下文"""
    logger.info(f"Building context from semantic results ({len(semantic_results)} results)...")
    
    context_parts = []
    for res in semantic_results:
        context_parts.append(f"- {res['text']} (Relevance: {res['similarity']:.2f})")
    
    context = "\n".join(context_parts)
    
    # 添加关键词结果作为补充
    if keyword_results and any(len(data) > 0 for data in keyword_results.values()):
        logger.info("Adding keyword results as supplement...")
        context += "\n\nAdditional database records:"
        for table_name, records in keyword_results.items():
            for rec in records[:2]:
                valid_fields = [
                    f"{k}: {v}" 
                    for k, v in rec.items() 
                    if v != "Not available" and "Sequence data is hidden" not in str(v)
                ]
                if valid_fields:
                    context += f"\n- {' | '.join(valid_fields)}"
    
    logger.info(f"Semantic context built, length: {len(context)} characters")
    logger.debug(f"Semantic context preview: {context[:200]}...")
    
    return context

def extract_entities(query: str):
    """提取实体（备用方案）"""
    logger.info(f"Extracting entities from query: '{query}'")
    
    pattern = r'\b[a-zA-Z0-9\-]{3,}\b'
    candidates = re.findall(pattern, query.lower())
    logger.debug(f"Entity candidates: {candidates}")
    
    stop_words = {
        'which', 'what', 'how', 'why', 'where', 'when', 'is', 'are', 'was', 'were',
        'does', 'do', 'did', 'the', 'and', 'or', 'for', 'to', 'of', 'in', 'on', 'with',
        'by', 'an', 'a', 'animal', 'hormones', 'gene', 'genes', 'related', 'receptor',
        'receptors', 'associated', 'function', 'role', 'effect', 'system'
    }
    
    filtered_terms = [term for term in candidates if term not in stop_words]
    logger.debug(f"Filtered terms: {filtered_terms}")
    
    if not filtered_terms:
        logger.info("No valid entities extracted")
        return None
    
    biological_terms = [t for t in filtered_terms if '-' in t or 'r' in t or 'l' in t or 'i' in t]
    entity = biological_terms[0] if biological_terms else filtered_terms[0]
    logger.info(f"Extracted entity: {entity}")
    
    return entity

# ==================== 主RAG函数（已修改以满足您的需求）====================
def rag_search(query: str):
    """主RAG搜索函数 - 本地数据库优先，始终补充网络搜索"""
    logger.info("=" * 50)
    logger.info(f"RAG Search started for query: '{query}'")
    start_time = datetime.now()
    
    # 1. 本地关键词检索
    logger.info("Step 1: Performing keyword retrieval...")
    keyword_results = keyword_retrieval(query)
    has_keyword_data = any(len(data) > 0 for data in keyword_results.values()) if keyword_results else False
    
    # 2. 判断查询类型并获取本地上下文
    is_simple = is_simple_query(query)
    logger.info(f"Query type: {'Simple' if is_simple else 'Complex'}")
    
    local_context = None
    local_context_type = "unknown"
    
    if is_simple and has_keyword_data:
        logger.info("Simple query with keyword data, using structured context")
        local_context = build_context_from_keyword(keyword_results)
        local_context_type = "structured"
    else:
        # 复杂查询或无关键词结果 - 尝试语义检索
        logger.info("Step 2: Performing semantic retrieval...")
        try:
            semantic_results = vector_store.semantic_search(query, top_k=3)
            logger.info(f"Semantic search returned {len(semantic_results)} results")
            
            if semantic_results:
                logger.info("Building context from semantic results")
                local_context = build_context_from_semantic(semantic_results, keyword_results)
                local_context_type = "semantic"
            elif has_keyword_data:
                logger.info("Semantic search returned no results, falling back to keyword results")
                local_context = build_context_from_keyword(keyword_results)
                local_context_type = "keyword_fallback"
            else:
                logger.info("No results from keyword or semantic search, trying entity extraction")
                entity = extract_entities(query)
                
                if entity:
                    logger.info(f"Retrying keyword search with entity: {entity}")
                    retry_results = keyword_retrieval(entity)
                    if any(len(data) > 0 for data in retry_results.values()):
                        local_context = build_context_from_keyword(retry_results)
                        local_context_type = "entity_extracted"
                        logger.info(f"Entity search successful, found context")
                    else:
                        local_context = f"No information found for '{entity}' in the database."
                        local_context_type = "no_results"
                        logger.info(f"Entity search returned no results")
                else:
                    local_context = "Could not find any relevant information in the database."
                    local_context_type = "no_results"
                    logger.info("No entity could be extracted")
                    
        except Exception as e:
            logger.error(f"Semantic search failed: {str(e)}")
            logger.error(traceback.format_exc())
            if has_keyword_data:
                logger.info("Falling back to keyword results due to semantic search error")
                local_context = build_context_from_keyword(keyword_results)
                local_context_type = "keyword_only"
            else:
                local_context = "Search failed. Please try a different query."
                local_context_type = "error"
    
    # 3. 无论本地是否有数据，都尝试进行网络搜索（补充信息）
    logger.info("Step 3: Performing web search for enrichment...")
    web_context = serper_search(query)
    if web_context:
        logger.info(f"Web search successful, retrieved context length: {len(web_context)}")
    else:
        logger.info("Web search returned no results or failed")
    
    # 4. 合并本地上下文和网络上下文
    final_context = None
    final_context_type = None
    
    if local_context and local_context_type not in ["no_results", "error"] and web_context:
        # 既有本地数据又有网络数据，合并
        logger.info("Merging local and web contexts")
        final_context = f"【Database Information】\n{local_context}\n\n【Web Search Information】\n{web_context}"
        final_context_type = "hybrid"
    elif local_context and local_context_type not in ["no_results", "error"]:
        # 只有本地数据
        logger.info("Using only local context")
        final_context = local_context
        final_context_type = local_context_type
    elif web_context:
        # 只有网络数据
        logger.info("Using only web context")
        final_context = web_context
        final_context_type = "web_search"
    else:
        # 都无数据
        logger.info("No context available from either source")
        final_context = local_context or "No relevant information found in the database or web search."
        final_context_type = local_context_type if local_context_type not in ["unknown"] else "no_results"
    
    # 5. 使用DeepSeek生成回答
    logger.info(f"Step 4: Generating answer with DeepSeek (final_context_type: {final_context_type})")
    logger.info(f"Final context length: {len(final_context) if final_context else 0} characters")
    
    if final_context and final_context_type not in ["no_results", "error"]:
        llm_result = deepseek_rag_client.generate_with_context(query, final_context, final_context_type)
        
        if llm_result["success"]:
            logger.info("DeepSeek generation successful")
            result = {
                "success": True,
                "answer": llm_result["response"],
                "model": llm_result["model"],
                "context_type": final_context_type,
                "has_data": True,
                "token_usage": llm_result.get("usage", {})
            }
        else:
            logger.error(f"DeepSeek generation failed: {llm_result['error']}")
            result = {
                "success": True,
                "answer": f"[Retrieved information but generation failed: {llm_result['error']}]\n\nRelevant information:\n{final_context}",
                "model": "rag-fallback",
                "context_type": final_context_type,
                "has_data": True,
                "generation_error": llm_result["error"]
            }
    else:
        logger.warning("No valid context to generate answer")
        result = {
            "success": True,
            "answer": final_context or "No relevant information found.",
            "model": "rag-no-data",
            "context_type": final_context_type,
            "has_data": False
        }
    
    # 计算耗时
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    logger.info(f"RAG search completed in {duration:.2f} seconds")
    logger.info(f"Result: has_data={result['has_data']}, context_type={result['context_type']}")
    logger.info(f"Answer preview: {result['answer'][:200]}...")
    logger.info("=" * 50)
    
    return result

# ==================== 工具函数 ====================
def check_system_status():
    """检查系统状态"""
    logger.info("Checking system status...")
    
    try:
        vector_count = vector_store.collection.count() if vector_store.collection else 0
        status = {
            "vector_store": {
                "initialized": True,
                "document_count": vector_count,
                "cache_size": len(vector_store.query_cache)
            },
            "deepseek": {
                "configured": True,
                "model": DEEPSEEK_MODEL,
                "api_key_valid": bool(DEEPSEEK_API_KEY and DEEPSEEK_API_KEY.startswith("sk-"))
            }
        }
        logger.info(f"System status: vector_count={vector_count}, cache_size={len(vector_store.query_cache)}")
        return status
    except Exception as e:
        logger.error(f"Error checking system status: {str(e)}")
        return {
            "vector_store": {"initialized": False, "error": str(e)},
            "deepseek": {"configured": True, "model": DEEPSEEK_MODEL}
        }

def force_rebuild_index():
    """强制重建索引"""
    logger.info("=" * 50)
    logger.info("Force rebuild index requested")
    
    try:
        logger.info(f"Deleting existing collection: {vector_store.collection.name}")
        collection_name = vector_store.collection.name
        vector_store.client.delete_collection(collection_name)
        
        logger.info("Creating new collection...")
        vector_store.collection = vector_store.client.create_collection(
            name=collection_name,
            metadata={"hnsw:space": "cosine"}
        )
        
        logger.info("Starting rebuild...")
        success = vector_store.rebuild_index()
        
        if success:
            logger.info("Index rebuild completed successfully")
            result = {"status": "success", "message": "Index rebuilt successfully"}
        else:
            logger.error("Index rebuild failed")
            result = {"status": "failed", "message": "Index rebuild failed"}
        
        logger.info("=" * 50)
        return result
        
    except Exception as e:
        logger.error(f"Error rebuilding index: {str(e)}")
        logger.error(traceback.format_exc())
        return {"status": "error", "message": str(e)}

# 启动时初始化
logger.info("Starting async vector store initialization...")
initialize_vector_store_async()