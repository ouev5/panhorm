from django.db.models import Q
from hormone_app.models import HormoneRelatedGene, HormoneReceptorFull, HormoneReceptorInfo
import numpy as np
from sentence_transformers import SentenceTransformer
import faiss
import threading
import re
from itertools import chain
# ==================== Keyword Retrieval Module ====================
def keyword_retrieval(query: str):
    """数据库关键词检索：修复表匹配+隐藏基因序列"""
    if not query or query.strip() == "":
        return None

    # 字段名映射：数据库字段 → 友好展示名
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
    
    # 需要隐藏的序列字段列表
    SEQUENCE_FIELDS = [
        "gene_sequence", 
        "protein_sequence", 
        "dna_sequence", 
        "mrna_sequence", 
        "amino_acid_sequence",
        "sequence",
        "receptor_coding_genes_sequence",
        "hormone_coding_genes_sequence"
    ]

    def format_item(item):
        """格式化条目：隐藏所有序列+添加提示文本"""
        formatted = {}
        for k, v in item.items():
            field_name = FIELD_MAPPING.get(k, k)
            # 隐藏所有序列字段，替换为提示文本
            if any(seq_field in k.lower() for seq_field in SEQUENCE_FIELDS) or "sequence" in k.lower():
                if v and isinstance(v, str) and len(v) > 0:
                    formatted[field_name] = "Sequence data is hidden for security"
                else:
                    formatted[field_name] = "Not available"
            else:
                # 其他字段完整显示（空值替换为Not available）
                formatted[field_name] = v if v is not None else "Not available"
        return formatted

    # 1. 检索HormoneRelatedGene表 (此表无receptor_name字段，保持原样)
    hgd_query = HormoneRelatedGene.objects.filter(
        Q(hormone_name__iexact=query) |
        Q(hormone_name__icontains=query) |
        Q(related_genes__icontains=query)
    )
    hgd_data = [format_item(item) for item in hgd_query.values()] if hgd_query.exists() else []

    # 2. 检索HormoneReceptorFull表 (新增对 receptor_name 的查询)
    phr_query = HormoneReceptorFull.objects.filter(
        Q(hormone_name__iexact=query) |
        Q(hormone_name__icontains=query) |
        Q(receptor_name__iexact=query) |  # 👈 新增
        Q(receptor_name__icontains=query) # 👈 新增
    )
    phr_data = [format_item(item) for item in phr_query.values()] if phr_query.exists() else []

    # 3. 检索HormoneReceptorInfo表 (新增对 receptor_name 的查询)
    nphr_query = HormoneReceptorInfo.objects.filter(
        Q(hormone_name__iexact=query) |
        Q(hormone_name__icontains=query) |
        Q(receptor_name__iexact=query) |  # 👈 新增
        Q(receptor_name__icontains=query) # 👈 新增
    )
    nphr_data = [format_item(item) for item in nphr_query.values()] if nphr_query.exists() else []

    # 强制返回所有表（含空表）
    final_result = {
        "hormone_related_gene": hgd_data,
        "hormone_receptor_full": phr_data,
        "hormone_receptor_info": nphr_data
    }
    
    return final_result

# ==================== 语义检索模块（启用） ====================
class HormoneVectorStore:
    def __init__(self, model_name="/www/wwwroot/default/animal_hormone/local_models/all-MiniLM-L6-v2/"):
        """初始化向量存储"""
        self.model_name = model_name
        self.model = None
        self.index = None
        self.texts = []  # 存储原始文本
        self.metadata = []  # 存储元数据（表名、ID等）
        self.is_loaded = False
        self.loading = False  # 防止重复加载
    
    def initialize_model(self):
        """延迟加载模型（节省启动时间）- 修复本地加载Pooling参数问题"""
        if self.model is None:
            try:
                print(f"正在加载句子嵌入模型: {self.model_name}")
                # 手动构建模型（解决Pooling参数缺失问题）
                from sentence_transformers import SentenceTransformer, models
                # 1. 加载预训练模型（本地路径）
                word_embedding_model = models.Transformer(self.model_name)
                # 2. 添加池化层（手动指定维度，解决参数缺失）
                pooling_model = models.Pooling(
                    word_embedding_model.get_word_embedding_dimension(),
                    pooling_mode_mean_tokens=True,
                    pooling_mode_cls_token=False,
                    pooling_mode_max_tokens=False
                )
                # 3. 组合模型
                self.model = SentenceTransformer(modules=[word_embedding_model, pooling_model])
                print("句子嵌入模型加载完成（手动构建，解决Pooling参数问题）")
            except Exception as e:
                print(f"加载句子嵌入模型失败: {e}")
                return False
        return True
    
    def load_data_to_vector(self):
        """从数据库加载数据并创建向量索引"""
        if self.loading or self.is_loaded:
            return
        
        self.loading = True
        print("开始构建语义检索向量索引...")
        
        if not self.initialize_model():
            self.loading = False
            return
        
        try:
            # 尝试加载预计算的向量索引（GPU服务器预计算）
            import os
            precompute_dir = "/www/wwwroot/default/animal_hormone"
            emb_path = os.path.join(precompute_dir, "rag_embeddings.npz")
            texts_path = os.path.join(precompute_dir, "rag_texts.json")
            meta_path = os.path.join(precompute_dir, "rag_metadata.json")
            
            if os.path.exists(emb_path) and os.path.exists(texts_path):
                print("检测到预计算向量索引，直接加载...")
                import json
                if not self.initialize_model():
                    self.loading = False
                    return
                
                data_emb = np.load(emb_path)
                embeddings = data_emb["embeddings"]
                
                with open(texts_path) as f:
                    all_texts = json.load(f)
                with open(meta_path) as f:
                    all_metadata = json.load(f)
                
                print(f"预计算数据: {len(all_texts)} 条文本, embeddings shape: {embeddings.shape}")
                
                # L2归一化（确保查询和数据库向量在同一尺度）
                norms = np.linalg.norm(embeddings, axis=1, keepdims=True)
                norms = np.maximum(norms, 1e-8)
                embeddings = embeddings / norms
                print(f"Embeddings L2 normalized, new norm: {np.linalg.norm(embeddings[0]):.4f}")
                
                # 直接创建FAISS索引
                dimension = embeddings.shape[1]
                self.index = faiss.IndexFlatL2(dimension)
                self.index.add(embeddings.astype("float32"))
                self.texts = all_texts
                self.metadata = all_metadata
                self.is_loaded = True
                print(f"预计算向量索引加载完成，包含 {len(all_texts)} 条记录")
                return
            
            # 没有预计算文件，走原始计算流程
            all_texts = []
            all_metadata = []
            
            # 1. 从HormoneRelatedGene表收集数据（含新字段）
            print("收集 HormoneRelatedGene 表数据...")
            for item in HormoneRelatedGene.objects.filter(organism__isnull=False).exclude(organism=""):
                if item.hormone_name:
                    org = getattr(item, 'organism', None) or ''
                    reg = getattr(item, 'regulation_type', None) or ''
                    
                    # 基础信息 + 物种
                    if item.related_genes and item.related_genes != "Not available":
                        text = f"Hormone {item.hormone_name} ({org}) {reg}regulates gene {item.related_genes}"
                        if item.pmid:
                            text += f" (PMID: {item.pmid})"
                        all_texts.append(text)
                        all_metadata.append({
                            'source': 'HormoneRelatedGene',
                            'id': item.id,
                            'field': 'related_genes'
                        })
                    
                    # 疾病关联 + DO标准化
                    if item.related_diseases and item.related_diseases != "Not available":
                        text = f"Gene {item.related_genes} is associated with disease: {item.related_diseases}"
                        do_name = getattr(item, 'doid_standardized_name', None) or ''
                        do_id = item.DO_ID or ''
                        if do_name:
                            text += f" (DO: {do_name}, {do_id})"
                        pmid_gd = getattr(item, 'pmid_gene_disease', None) or ''
                        if pmid_gd:
                            text += f" (PMID: {pmid_gd})"
                        all_texts.append(text)
                        all_metadata.append({
                            'source': 'HormoneRelatedGene',
                            'id': item.id,
                            'field': 'related_diseases'
                        })
            
            # 2. 从HormoneReceptorFull表收集数据
            print("收集 HormoneReceptorFull 表数据...")
            for item in HormoneReceptorFull.objects.all():
                if item.hormone_name and item.receptor_name:
                    text = f"Hormone {item.hormone_name} binds to receptor {item.receptor_name}"
                    all_texts.append(text)
                    all_metadata.append({
                        'source': 'HormoneReceptorFull',
                        'id': item.id,
                        'field': 'receptor_binding'
                    })
            
            # 3. 从HormoneReceptorInfo表收集数据
            print("收集 HormoneReceptorInfo 表数据...")
            for item in HormoneReceptorInfo.objects.all():
                if item.hormone_name and item.receptor_name:
                    text = f"Hormone {item.hormone_name} interacts with receptor {item.receptor_name}"
                    if item.pubchem_id and item.pubchem_id != "Not available":
                        text += f" (PubChem ID: {item.pubchem_id})"
                    all_texts.append(text)
                    all_metadata.append({
                        'source': 'HormoneReceptorInfo',
                        'id': item.id,
                        'field': 'receptor_interaction'
                    })
            
            if not all_texts:
                print("警告：没有数据可用于构建向量索引")
                self.loading = False
                return
            
            print(f"收集到 {len(all_texts)} 条文本用于向量化...")
            
            # 生成向量嵌入（分批处理避免内存问题）
            batch_size = 64
            embeddings_list = []
            
            for i in range(0, len(all_texts), batch_size):
                batch_texts = all_texts[i:i+batch_size]
                try:
                    batch_embeddings = self.model.encode(batch_texts, show_progress_bar=False)
                    embeddings_list.append(batch_embeddings)
                    processed = min(i + batch_size, len(all_texts))
                    print(f"已处理 {processed}/{len(all_texts)} 条文本")
                except Exception as e:
                    print(f"处理批次 {i//batch_size + 1} 时出错: {e}")
                    continue
            
            if not embeddings_list:
                print("错误：未能生成任何向量嵌入")
                self.loading = False
                return
            
            # 合并所有嵌入
            embeddings = np.vstack(embeddings_list)
            
            # 创建FAISS索引
            dimension = embeddings.shape[1]
            self.index = faiss.IndexFlatL2(dimension)
            self.index.add(embeddings.astype('float32'))
            
            self.texts = all_texts
            self.metadata = all_metadata
            self.is_loaded = True
            
            print(f"语义检索向量索引构建完成，包含 {len(all_texts)} 条记录")
            
        except Exception as e:
            print(f"构建向量索引时出错: {e}")
        finally:
            self.loading = False
    
    def semantic_search(self, query: str, top_k: int = 30):
        """语义搜索：返回最相关的top_k个结果"""
        if not self.is_loaded or not self.index or not self.model:
            return []
        
        try:
            # 编码查询
            query_embedding = self.model.encode([query], show_progress_bar=False)
            # L2归一化查询向量
            import numpy as _np
            _norm = _np.linalg.norm(query_embedding, axis=1, keepdims=True)
            _norm = _np.maximum(_norm, 1e-8)
            query_embedding = query_embedding / _norm
            
            # 搜索最相似的k个（多检索一些用于去重）
            search_k = min(top_k * 3, len(self.texts))
            distances, indices = self.index.search(query_embedding.astype('float32'), search_k)
            
            # 返回结果（去重相似的文本）
            results = []
            seen_sources = set()
            
            for idx, distance in zip(indices[0], distances[0]):
                if idx < len(self.texts):
                    metadata = self.metadata[idx]
                    source_key = f"{metadata['source']}_{metadata['id']}"
                    
                    # 去重：同一个数据库记录只返回一次
                    if source_key not in seen_sources:
                        seen_sources.add(source_key)
                        
                        # 转换为相似度分数（0-1）
                        similarity = float(1 / (1 + distance))
                        
                        results.append({
                            'text': self.texts[idx],
                            'similarity': similarity,
                            'metadata': metadata
                        })
                        
                        if len(results) >= top_k:
                            break
            
            return results
            
        except Exception as e:
            print(f"语义搜索时出错: {e}")
            return []


# ==================== 初始化向量存储 ====================
# 创建全局向量存储实例
vector_store = HormoneVectorStore()

def initialize_vector_store_async():
    """异步初始化向量存储（在后台线程中运行）"""
    def load_in_background():
        try:
            vector_store.load_data_to_vector()
        except Exception as e:
            print(f"后台加载向量索引失败: {e}")
    
    # 启动后台线程加载向量索引
    thread = threading.Thread(target=load_in_background, daemon=True)
    thread.start()
    print("已启动后台线程加载语义检索索引...")


# ==================== 混合检索+Prompt构建 ====================
def is_short_query(query: str) -> bool:
    """判断是否为短词查询（基因/激素名）"""
    query = query.strip()
    if len(query) == 0:
        return False
    # 短于4个词，且不含问号、疑问词
    words = query.split()
    if len(words) > 3:
        return False
    if any(w.lower() in ['which', 'what', 'how', 'why', 'where', 'is', 'are', 'does'] for w in words):
        return False
    if '?' in query or '.' in query:
        return False
    return True


def hybrid_advanced_retrieval(query: str):
    """混合检索：支持短词结构化返回 + 长句语义生成"""
    short_query = is_short_query(query)
    
    # 1. 关键词检索（始终执行）
    db_data = keyword_retrieval(query)
    has_keyword_data = any(len(data) > 0 for data in db_data.values()) if db_data else False
    
    # 2. 语义检索（仅长句时启用）
    semantic_results = []
    semantic_status = "disabled"
    
    if not short_query and vector_store.is_loaded:
        semantic_results = vector_store.semantic_search(query, top_k=30)
        semantic_status = "enabled"
    elif vector_store.loading:
        semantic_status = "loading"
    
    return {
        "has_db_data": has_keyword_data or len(semantic_results) > 0,
        "keyword_content": db_data if db_data else {},
        "semantic_content": semantic_results,
        "is_short_query": short_query,  # 👈 新增标志
        "semantic_status": semantic_status
    }


def build_biogpt_prompt(query: str, retrieval_result: dict):
    is_short = retrieval_result.get("is_short_query", False)
    
    if is_short:
        # ========== 短词：直接返回关键词检索结果 ==========
        lines = []
        keyword_data = retrieval_result["keyword_content"]
        for table_name, records in keyword_data.items():
            for rec in records:
                valid_fields = [
                    f"{k}: {v}" 
                    for k, v in rec.items() 
                    if v != "Not available" and "Sequence data is hidden" not in str(v)
                ]
                if valid_fields:
                    lines.append(" | ".join(valid_fields))
        content = "\n".join(lines) if lines else "No matching hormone, gene, or receptor found in the database."
        return {"type": "structured", "content": content}
    
    else:
        # ========== 长句：优先使用语义检索结果 ==========
        semantic_results = retrieval_result.get("semantic_content", [])
        
        if semantic_results:
            extracted_pairs = set()
            for res in semantic_results:
                text = res['text']
                # 从预定义的文本模板中提取 hormone 和 gene/receptor
                # 模板1: "Hormone X is associated with genes: Y"
                if "is associated with genes:" in text:
                    try:
                        parts = text.split("is associated with genes:")
                        hormone = parts[0].replace("Hormone ", "").strip()
                        genes = parts[1].strip().rstrip('.')
                        # 支持多个基因（用逗号分隔）
                        for gene in [g.strip() for g in genes.split(',')]:
                            if hormone and gene:
                                extracted_pairs.add(f"Hormone: {hormone} | Gene: {gene}")
                    except:
                        pass
                
                # 模板2: "Hormone X binds to receptor Y"
                elif "binds to receptor" in text:
                    try:
                        parts = text.split("binds to receptor")
                        hormone = parts[0].replace("Hormone ", "").strip()
                        receptor = parts[1].strip().rstrip('.')
                        if hormone and receptor:
                            extracted_pairs.add(f"Hormone: {hormone} | Gene: {receptor}")
                    except:
                        pass
                
                # 模板3: "Hormone X interacts with receptor Y"
                elif "interacts with receptor" in text:
                    try:
                        parts = text.split("interacts with receptor")
                        hormone = parts[0].replace("Hormone ", "").strip()
                        receptor = parts[1].split(" (PubChem")[0].strip().rstrip('.')
                        if hormone and receptor:
                            extracted_pairs.add(f"Hormone: {hormone} | Gene: {receptor}")
                    except:
                        pass
            
            if extracted_pairs:
                return {
                    "type": "structured",
                    "content": "\n".join(sorted(extracted_pairs))
                }
        
        # ========== 备用方案：如果语义检索无结果，回退到实体提取+精确匹配 ==========
        import re
        pattern = r'\b[a-zA-Z0-9\-]{3,}\b'
        candidates = re.findall(pattern, query.lower())
        
        stop_words = {
            'which', 'what', 'how', 'why', 'where', 'when', 'is', 'are', 'was', 'were',
            'does', 'do', 'did', 'the', 'and', 'or', 'for', 'to', 'of', 'in', 'on', 'with',
            'by', 'an', 'a', 'animal', 'hormones', 'gene', 'genes', 'related', 'receptor',
            'receptors', 'associated', 'function', 'role', 'effect', 'system'
        }
        filtered_terms = [term for term in candidates if term not in stop_words]
        if not filtered_terms:
            return {"type": "structured", "content": "Could not extract a valid gene or hormone name from the question."}
        
        biological_terms = [t for t in filtered_terms if '-' in t or 'r' in t or 'l' in t or 'i' in t]
        target_gene = biological_terms[0] if biological_terms else filtered_terms[0]
        target_normalized = target_gene.replace('-', '').replace(' ', '').lower()
        
        all_hormone_gene_pairs = set()
        for item in HormoneRelatedGene.objects.filter(organism__isnull=False).exclude(organism=""):
            if item.hormone_name and item.related_genes:
                db_norm = str(item.related_genes).replace('-', '').replace(' ', '').lower()
                if db_norm == target_normalized:
                    all_hormone_gene_pairs.add(f"Hormone: {item.hormone_name} | Gene: {item.related_genes}")
        
        for item in chain(HormoneReceptorFull.objects.all(), HormoneReceptorInfo.objects.all()):
            if item.hormone_name and item.receptor_name:
                db_norm = str(item.receptor_name).replace('-', '').replace(' ', '').lower()
                if db_norm == target_normalized:
                    all_hormone_gene_pairs.add(f"Hormone: {item.hormone_name} | Gene: {item.receptor_name}")
        
        if all_hormone_gene_pairs:
            content = "\n".join(sorted(all_hormone_gene_pairs))
        else:
            content = f"No animal hormones found that are directly associated with the gene/receptor '{target_gene}'."
        
        return {"type": "structured", "content": content}
# ==================== 工具函数 ====================
def check_vector_store_status():
    """检查向量存储状态"""
    return {
        "is_loaded": vector_store.is_loaded,
        "is_loading": vector_store.loading,
        "model_name": vector_store.model_name if hasattr(vector_store, 'model_name') else None,
        "text_count": len(vector_store.texts) if vector_store.texts else 0
    }


def force_reload_vector_store():
    """强制重新加载向量存储（用于数据库更新后）"""
    if vector_store.loading:
        return {"status": "already_loading", "message": "向量索引正在加载中"}
    
    # 重置状态
    vector_store.is_loaded = False
    vector_store.texts = []
    vector_store.metadata = []
    vector_store.index = None
    
    # 重新加载
    #initialize_vector_store_async()
    #return {"status": "reloading", "message": "已开始重新加载向量索引"}