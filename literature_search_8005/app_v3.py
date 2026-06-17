"""
文献检索服务 v3.0 - 异步任务模式
后台分析完保存到文件，前端轮询加载，数据全量保留
"""

import os, re, json, time, uuid, threading, requests
from typing import Optional, List, Dict
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import HTMLResponse, FileResponse
from starlette.middleware.gzip import GZipMiddleware
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
import uvicorn

DEEPSEEK_API_KEY = "YOUR_VALUE_HERE"
DEEPSEEK_API_URL = "https://api.deepseek.com/v1/chat/completions"
DEEPSEEK_MODEL = "deepseek-v4-pro"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_EMAIL = "bioinfo@example.com"
RESULTS_DIR = "/www/wwwroot/literature_search/results"
os.makedirs(RESULTS_DIR, exist_ok=True)
tasks = {}
def call_deepseek(system_prompt: str, user_prompt: str, temperature: float = 0.3, max_tokens: int = 4000, timeout: int = 180) -> str:
    """通用DeepSeek调用"""
    headers = {"Content-Type": "application/json", "Authorization": f"Bearer {DEEPSEEK_API_KEY}"}
    data = {
        "model": DEEPSEEK_MODEL,
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt}
        ],
        "temperature": temperature,
        "max_tokens": max_tokens,
        "thinking": {"type": "disabled"}
    }
    resp = requests.post(DEEPSEEK_API_URL, headers=headers, json=data, timeout=timeout)
    resp.raise_for_status()
    return resp.json()["choices"][0]["message"]["content"]


def optimize_query(user_query: str) -> dict:
    """让DeepSeek优化用户的检索词，生成PubMed最优检索式"""
    system = """You are a PubMed search expert. Translate user query to English and create optimal PubMed search terms.

CRITICAL: Return ONLY a valid JSON object. Do NOT use double quotes inside any string value - use single quotes or no quotes instead.
Example: {"translated": "insulin related genes", "search_term": "insulin[MeSH] AND genes[MeSH]", "related_terms": ["insulin signaling genes", "insulin receptor gene"]}

Fields:
- translated: English translation
- search_term: PubMed Boolean query using AND/OR, field tags like [MeSH], [Title/Abstract], but NO double quotes inside the value
- related_terms: array of 4 related search queries

Return ONLY the JSON object, nothing else."""

    user = f'User query: "{user_query}"'
    result = call_deepseek(system, user, temperature=0.2, max_tokens=1000)
    
    # 灵活解析JSON
    try:
        cleaned = re.sub(r'```(?:json)?\s*', '', result).strip()
        json_match = re.search(r'\{[\s\S]*\}', cleaned)
        if json_match:
            json_str = json_match.group()
            # 先尝试标准解析
            try:
                data = json.loads(json_str)
            except json.JSONDecodeError:
                # 用ast.literal_eval作为备选
                import ast
                try:
                    data = ast.literal_eval(json_str)
                except:
                    # 最后：尝试修复引号问题后重新解析
                    # 把检索式内部的"替换为'
                    fixed = re.sub(r'(?<=[a-zA-Z\[])\"(?=[a-zA-Z\]])', "'", json_str)
                    try:
                        data = json.loads(fixed)
                    except:
                        data = {}
            
            # 兼容不同键名
            translated = data.get("translated") or data.get("translation") or data.get("translated_text") or data.get("english") or user_query
            search_term = data.get("search_term") or data.get("search") or data.get("query") or data.get("pubmed_query") or user_query
            related = data.get("related_terms") or data.get("related") or data.get("alternatives") or []
            return {
                "translated": translated,
                "search_term": search_term,
                "related_terms": related
            }
    except Exception as e:
        print(f"optimize_query parse error: {e}, raw: {result[:200]}")
    
    # fallback
    return {
        "translated": user_query,
        "search_term": user_query,
        "related_terms": []
    }


def deep_search_queries(user_query: str) -> list:
    """深度检索：让DeepSeek生成多个不同角度的检索词"""
    system = """You are a PubMed search expert. Generate multiple diverse search queries to comprehensively cover a research topic.
You MUST return ONLY valid JSON (no markdown, no code blocks, no extra text).
The JSON must be an object with a "queries" array. Each element has "label" and "search_term".
IMPORTANT: In search_term values, do NOT use double quotes inside the string. Use single quotes or no quotes instead.
Example format: {"queries": [{"label": "Core query", "search_term": "PCOS AND biomarkers AND diagnosis"}, {"label": "Mechanism", "search_term": "PCOS AND insulin resistance AND pathway"}]}

Generate exactly 5 queries from DIFFERENT angles:
1. Core/direct query
2. Mechanism/pathway angle
3. Disease/phenotype angle
4. Comparative/evolutionary angle
5. Clinical/therapeutic angle

Use PubMed Boolean operators (AND, OR) but NO double quotes inside search_term values.
Return ONLY the JSON object."""

    user = f'User query: "{user_query}"\n\nGenerate 5 diverse PubMed search queries. Return ONLY JSON.'
    result = call_deepseek(system, user, temperature=0.3, max_tokens=2000)
    try:
        cleaned = re.sub(r'```(?:json)?\s*', '', result).strip()
        # 如果返回的是数组格式，包装成对象
        if cleaned.startswith('['):
            cleaned = '{"queries": ' + cleaned + '}'
        json_match = re.search(r'\{[\s\S]*\}', cleaned)
        if json_match:
            json_str = json_match.group()
            try:
                data = json.loads(json_str)
            except json.JSONDecodeError:
                import ast
                try:
                    data = ast.literal_eval(json_str)
                except:
                    # 最后尝试：用正则提取每个query
                    queries = []
                    label_matches = re.finditer(r'"label"\s*:\s*"([^"]*)"', json_str)
                    term_matches = re.finditer(r'"search_term"\s*:\s*"((?:[^"\\]|\\.)*)"', json_str)
                    labels = [m.group(1) for m in label_matches]
                    terms = [m.group(1).replace('\\"', '"') for m in term_matches]
                    for i in range(min(len(labels), len(terms))):
                        queries.append({"label": labels[i], "search_term": terms[i]})
                    if queries:
                        return queries
                    print(f"deep_search_queries all parse methods failed, raw: {json_str[:200]}")
                    return []
            queries = data.get("queries", [])
            # 兼容：如果元素里没有search_term但有query字段
            for q in queries:
                if "search_term" not in q and "query" in q:
                    q["search_term"] = q["query"]
                if "label" not in q and "description" in q:
                    q["label"] = q["description"][:30]
            return queries
    except Exception as e:
        print(f"deep_search_queries parse error: {e}")
    return []


# ==================== PubMed检索模块 ====================

class PubMedSearcher:
    """PubMed文献检索器"""
    
    def search(self, query: str, max_results: int = 15) -> List[Dict]:
        """搜索PubMed，返回文献列表"""
        try:
            params = {
                "db": "pubmed",
                "term": query,
                "retmax": max_results,
                "sort": "relevance",
                "retmode": "json",
                "email": NCBI_EMAIL
            }
            resp = requests.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            pmid_list = data.get("esearchresult", {}).get("idlist", [])
            
            if not pmid_list:
                return []
            
            fetch_params = {
                "db": "pubmed",
                "id": ",".join(pmid_list),
                "rettype": "xml",
                "retmode": "xml",
                "email": NCBI_EMAIL
            }
            fetch_resp = requests.get(f"{EUTILS_BASE}/efetch.fcgi", params=fetch_params, timeout=30)
            fetch_resp.raise_for_status()
            
            articles = self._parse_xml(fetch_resp.text, pmid_list)
            return articles
            
        except Exception as e:
            print(f"PubMed搜索失败: {e}")
            return []
    
    def _parse_xml(self, xml_text: str, pmid_list: List[str]) -> List[Dict]:
        """解析PubMed XML响应"""
        articles = []
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError:
            return []
        
        for article_elem in root.findall(".//PubmedArticle"):
            try:
                medline = article_elem.find("MedlineCitation")
                if medline is None:
                    continue
                article_data = medline.find("Article")
                if article_data is None:
                    continue
                
                title_elem = article_data.find("ArticleTitle")
                title = title_elem.text if title_elem is not None and title_elem.text else "N/A"
                
                author_list = article_data.find("AuthorList")
                authors = []
                if author_list is not None:
                    for author in author_list.findall("Author")[:5]:
                        last = author.find("LastName")
                        initials = author.find("Initials")
                        last_text = last.text if last is not None and last.text else ""
                        init_text = initials.text if initials is not None and initials.text else ""
                        if last_text:
                            authors.append(f"{last_text} {init_text}")
                author_str = ", ".join(authors)
                if author_list is not None and len(author_list.findall("Author")) > 5:
                    author_str += " et al."
                
                journal_elem = article_data.find("Journal")
                journal = ""
                journal_abbrev = ""
                year = ""
                if journal_elem is not None:
                    title_j = journal_elem.find("Title")
                    iso = journal_elem.find("ISOAbbreviation")
                    journal = title_j.text if title_j is not None and title_j.text else ""
                    journal_abbrev = iso.text if iso is not None and iso.text else journal
                    pub_date = journal_elem.find("PubDate")
                    if pub_date is not None:
                        year_elem = pub_date.find("Year")
                        if year_elem is not None and year_elem.text:
                            year = year_elem.text
                
                abstract_parts = []
                abstract_elem = article_data.find("Abstract")
                if abstract_elem is not None:
                    for text in abstract_elem.findall("AbstractText"):
                        if text.text:
                            label = text.get("Label", "")
                            if label:
                                abstract_parts.append(f"{label}: {text.text}")
                            else:
                                abstract_parts.append(text.text)
                abstract = " ".join(abstract_parts)[:500]
                
                pmid_elem = medline.find("PMID")
                pmid = pmid_elem.text if pmid_elem is not None and pmid_elem.text else ""
                
                doi = ""
                for eid in article_data.findall("ELocationID"):
                    if eid.get("EIdType") == "doi" and eid.text:
                        doi = eid.text
                        break
                
                mesh_terms = []
                mesh_list = medline.find("MeshHeadingList")
                if mesh_list is not None:
                    for mesh in mesh_list.findall("MeshHeading")[:10]:
                        desc = mesh.find("DescriptorName")
                        if desc is not None and desc.text:
                            mesh_terms.append(desc.text)
                
                articles.append({
                    "pmid": pmid,
                    "title": title,
                    "authors": author_str,
                    "journal": journal_abbrev,
                    "year": year,
                    "abstract": abstract,
                    "doi": doi,
                    "mesh_terms": mesh_terms,
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                })
            except Exception as e:
                print(f"解析文章失败: {e}")
                continue
        return articles
    
    def search_genes(self, query: str, max_results: int = 20) -> List[Dict]:
        """搜索基因数据库"""
        try:
            params = {
                "db": "gene",
                "term": f"{query}[Title] AND Homo sapiens[Organism]",
                "retmax": max_results,
                "retmode": "json",
                "email": NCBI_EMAIL
            }
            resp = requests.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            gene_ids = data.get("esearchresult", {}).get("idlist", [])
            if not gene_ids:
                return []
            sum_params = {
                "db": "gene",
                "id": ",".join(gene_ids[:10]),
                "retmode": "json",
                "email": NCBI_EMAIL
            }
            sum_resp = requests.get(f"{EUTILS_BASE}/esummary.fcgi", params=sum_params, timeout=30)
            sum_resp.raise_for_status()
            sum_data = sum_resp.json()
            genes = []
            result = sum_data.get("result", {})
            for gid in gene_ids[:10]:
                gene_info = result.get(gid, {})
                if gene_info and not isinstance(gene_info, str):
                    genes.append({
                        "gene_id": gid,
                        "symbol": gene_info.get("name", ""),
                        "description": gene_info.get("description", ""),
                        "url": f"https://www.ncbi.nlm.nih.gov/gene/{gid}"
                    })
            return genes
        except Exception as e:
            print(f"Gene搜索失败: {e}")
            return []


# ==================== AI分析模块 ====================

class AIAnalyzer:
    """AI分析器 - 并发版"""
    
    def _analyze_single_article(self, query: str, article: Dict) -> Dict:
        """分析单篇文章提取基因"""
        prompt = f"""Analyze this biomedical article about {query}.

Title: {article['title']}
Abstract: {article['abstract']}
MeSH: {', '.join(article.get('mesh_terms',[])[:5])}

TASK: Extract gene/protein symbols directly related to {query} from THIS article. Extract genes ONLY if:
1. The gene symbol is explicitly mentioned in the text (e.g. AMH, SHBG, INS, CYP17A1)
2. A protein/hormone is mentioned that maps to a specific gene (e.g. "anti-mullerian hormone" -> AMH, "insulin receptor" -> INSR)
3. A specific pathway is discussed with clear gene involvement (e.g. "insulin signaling pathway" when the article actually discusses insulin mechanism)
4. MeSH terms that directly name specific genes

DO NOT add generic pathway genes just because the topic is related. Only include genes with DIRECT textual evidence in this article.
Return max 5 most relevant genes per article. If no genes can be extracted, return empty genes array but still provide key_finding.

Respond with ONLY a JSON object (no other text):
{{"genes":[{{"symbol":"GENE1","full_name":"Full Name","role":"role description","finding":"what was found"}}],"key_finding":"main finding"}}"""
        try:
            result = call_deepseek(
                "You are an expert molecular biologist. Aggressively extract and INFER all gene/protein symbols from biomedical text. Map proteins/hormones to their encoding genes. Map pathways to key participant genes. Return ONLY valid JSON.",
                prompt, temperature=0.2, max_tokens=1500, timeout=60
            )
            result = result.strip()
            # 清理各种包裹
            for prefix in ["```json", "```", "```JSON"]:
                if result.startswith(prefix): result = result[len(prefix):]
            if result.endswith("```"): result = result[:-3]
            result = result.strip()
            # 找到第一个{和最后一个}
            start = result.find("{")
            end = result.rfind("}") + 1
            if start >= 0 and end > start:
                result = result[start:end]
            parsed = json.loads(result)
            for g in parsed.get("genes", []):
                g["pmid"] = article["pmid"]
                g["article_title"] = article["title"]
            return parsed
        except Exception as e:
            return {"genes": [], "key_finding": "", "_error": str(e)[:100]}
    
    def extract_genes_from_literature(self, query: str, articles: List[Dict], task_id: str = None) -> Dict:
        """并发分析所有文献 - 每篇文章独立调用AI"""
        all_genes = {}
        all_findings = []
        
        total = len(articles)
        with ThreadPoolExecutor(max_workers=3) as executor:
            futures = {executor.submit(self._analyze_single_article, query, art): art for art in articles}
            done_count = 0
            for future in as_completed(futures):
                done_count += 1
                if task_id and task_id in tasks:
                    tasks[task_id]["progress"] = 30 + int(done_count / total * 65)
                    tasks[task_id]["message"] = f"AI analyzing article {done_count}/{total}..."
                try:
                    result = future.result(timeout=90)
                    for g in result.get("genes", []):
                        symbol = g.get("symbol", "").upper()
                        if not symbol: continue
                        if symbol not in all_genes:
                            all_genes[symbol] = {"symbol": symbol, "full_name": g.get("full_name", ""), "role": g.get("role", ""), "evidence": []}
                        all_genes[symbol]["evidence"].append({"pmid": g.get("pmid", futures[future]["pmid"])})
                    if result.get("key_finding"):
                        all_findings.append(result["key_finding"])
                except: pass
        
        for sym in all_genes:
            sorted_genes = sorted(all_genes.values(), key=lambda x: len(x.get("evidence", [])), reverse=True)
        summary = f"Found {len(sorted_genes)} genes related to {query} across {len(articles)} articles."
        if sorted_genes:
            top = sorted_genes[0]
            summary += f" The most frequently mentioned gene is {top['symbol']} ({top.get('full_name','')})."
        
        return {"summary": summary, "genes": sorted_genes, "key_findings": all_findings[:5]}

# ==================== FastAPI服务 ====================



# ==================== 异步任务引擎 ====================

def _run_normal_search(task_id, query, max_results=15):
    try:
        tasks[task_id]["status"] = "running"
        tasks[task_id]["message"] = "AI optimizing search terms..."
        tasks[task_id]["progress"] = 5
        optimized = optimize_query(query)
        tasks[task_id]["message"] = "Searching PubMed..."
        tasks[task_id]["progress"] = 15
        articles = searcher.search(optimized.get("search_term", query), max_results)
        if not articles:
            _save_result(task_id, {"query":query,"optimized_query":optimized,"articles":[],"gene_analysis":None,"message":"未找到相关文献"})
            return
        tasks[task_id]["message"] = f"AI analyzing {len(articles)} articles..."
        tasks[task_id]["progress"] = 25
        ga = analyzer.extract_genes_from_literature(query, articles, task_id)
        _save_result(task_id, {"query":query,"optimized_query":optimized,"total_articles":len(articles),"articles":articles,"gene_analysis":ga})
    except Exception as e:
        tasks[task_id]["status"] = "failed"
        tasks[task_id]["message"] = str(e)[:200]

def _run_deep_search(task_id, query, max_per_query=10):
    try:
        tasks[task_id]["status"] = "running"
        tasks[task_id]["message"] = "AI generating multi-angle search terms..."
        tasks[task_id]["progress"] = 5
        queries = deep_search_queries(query)
        if not queries: queries = [{"label":"Default","search_term":query}]
        tasks[task_id]["message"] = f"Searching PubMed with {len(queries)} angles..."
        tasks[task_id]["progress"] = 15
        all_articles, seen, details = [], set(), []
        for i,q in enumerate(queries):
            term = q.get("search_term","")
            if not term: continue
            arts = searcher.search(term, max_per_query)
            nc = 0
            for a in arts:
                if a["pmid"] not in seen:
                    seen.add(a["pmid"])
                    all_articles.append(a)
                    nc += 1
            details.append({"label":q.get("label",""),"search_term":term,"found":len(arts),"new_unique":nc})
            tasks[task_id]["progress"] = 15 + int((i+1)/len(queries)*15)
        if not all_articles:
            _save_result(task_id, {"query":query,"search_details":details,"articles":[],"gene_analysis":None,"message":"未找到相关文献"})
            return
        tasks[task_id]["message"] = f"AI analyzing {len(all_articles)} articles for genes..."
        tasks[task_id]["progress"] = 30
        ga = analyzer.extract_genes_from_literature(query, all_articles, task_id)
        _save_result(task_id, {"query":query,"total_unique_articles":len(all_articles),"search_details":details,"articles":all_articles,"gene_analysis":ga})
    except Exception as e:
        tasks[task_id]["status"] = "failed"
        tasks[task_id]["message"] = str(e)[:200]

def _save_result(task_id, result):
    fp = os.path.join(RESULTS_DIR, f"{task_id}.json")
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
    tasks[task_id]["status"] = "done"
    tasks[task_id]["progress"] = 100
    tasks[task_id]["message"] = "Analysis complete!"
    tasks[task_id]["result_file"] = fp

app = FastAPI(title="文献检索服务 v3.0 (异步任务模式)")
app.add_middleware(GZipMiddleware, minimum_size=500)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_credentials=True, allow_methods=["*"], allow_headers=["*"])
searcher = PubMedSearcher()
analyzer = AIAnalyzer()

class SearchRequest(BaseModel):
    query: str
    max_results: int = 15

class DeepSearchRequest(BaseModel):
    query: str
    max_results_per_query: int = 10

@app.post("/api/search/submit")
def submit_search(request: SearchRequest):
    if not request.query.strip(): raise HTTPException(400, detail="查询不能为空")
    task_id = str(uuid.uuid4())[:8]
    tasks[task_id] = {"status":"queued","progress":0,"message":"Task created","start_time":time.time(),"mode":"normal"}
    threading.Thread(target=_run_normal_search, args=(task_id,request.query,request.max_results), daemon=True).start()
    return {"task_id":task_id,"status":"queued"}

@app.post("/api/deep-search/submit")
def submit_deep_search(request: DeepSearchRequest):
    if not request.query.strip(): raise HTTPException(400, detail="查询不能为空")
    task_id = str(uuid.uuid4())[:8]
    tasks[task_id] = {"status":"queued","progress":0,"message":"Task created","start_time":time.time(),"mode":"deep"}
    threading.Thread(target=_run_deep_search, args=(task_id,request.query,request.max_results_per_query), daemon=True).start()
    return {"task_id":task_id,"status":"queued"}

@app.get("/api/task/{task_id}")
def get_task_status(task_id: str):
    if task_id not in tasks: raise HTTPException(404, detail="Task not found")
    t = tasks[task_id]
    resp = {"task_id":task_id,"status":t["status"],"progress":t["progress"],"message":t["message"]}
    if t["status"] == "done": resp["result_url"] = f"/api/result/{task_id}"
    return resp

@app.get("/api/result/{task_id}")
def get_result(task_id: str):
    fp = os.path.join(RESULTS_DIR, f"{task_id}.json")
    if not os.path.exists(fp): raise HTTPException(404, detail="Result not found")
    return FileResponse(fp, media_type="application/json")

@app.get("/api/health")
def health():
    return {"status":"ok","service":"文献检索服务 v3.0 (异步)","active":len([t for t in tasks.values() if t["status"]=="running"]),"done":len([t for t in tasks.values() if t["status"]=="done"])}

@app.get("/")
async def root():
    return HTMLResponse(content=HTML_PAGE)
HTML_PAGE = r"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GeneLit - Literature Search & Gene Discovery</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap" rel="stylesheet">
    <link href="https://cdn.jsdelivr.net/npm/font-awesome@4.7.0/css/font-awesome.min.css" rel="stylesheet">
    <style>
        * { font-family: 'Inter', sans-serif; }
        .glass {
            background: #ffffff;
            border: 1px solid #e5e7eb;
            box-shadow: 0 1px 3px rgba(0,0,0,0.06);
        }
        .glass-hover:hover {
            box-shadow: 0 4px 12px rgba(139,92,246,0.12);
            border-color: rgba(139,92,246,0.3);
        }
        .pulse-glow { animation: pulseGlow 2s ease-in-out infinite; }
        @keyframes pulseGlow {
            0%, 100% { box-shadow: 0 0 15px rgba(139,92,246,0.1); }
            50% { box-shadow: 0 0 30px rgba(139,92,246,0.2); }
        }
        .gene-tag {
            background: linear-gradient(135deg, rgba(139,92,246,0.1), rgba(59,130,246,0.1));
            border: 1px solid rgba(139,92,246,0.25);
        }
        .evidence-line { border-left: 3px solid rgba(139,92,246,0.4); }
        .fade-in { animation: fadeIn 0.5s ease-out; }
        @keyframes fadeIn {
            from { opacity: 0; transform: translateY(10px); }
            to { opacity: 1; transform: translateY(0); }
        }
        .loader {
            border: 3px solid #e5e7eb;
            border-top: 3px solid #7c3aed;
            border-radius: 50%;
            width: 40px;
            height: 40px;
            animation: spin 1s linear infinite;
        }
        @keyframes spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
        }
        ::-webkit-scrollbar { width: 6px; }
        ::-webkit-scrollbar-track { background: #f9fafb; }
        ::-webkit-scrollbar-thumb { background: rgba(139,92,246,0.2); border-radius: 3px; }
        nav {
            background: linear-gradient(135deg, #0f172a, #1e3a8a);
            padding: 1rem 1.5rem;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            position: sticky; top: 0; z-index: 100;
        }
        .nav-container {
            max-width: 1400px; margin: 0 auto;
            display: flex; justify-content: space-between; align-items: center;
        }
        .back-btn {
            background: rgba(255,255,255,0.1); color: white;
            border: 1px solid rgba(255,255,255,0.2);
            padding: 0.6rem 1.2rem; border-radius: 8px;
            text-decoration: none; display: inline-flex;
            align-items: center; gap: 0.5rem; transition: all 0.2s;
        }
        .back-btn:hover { background: rgba(255,255,255,0.2); }
        .deep-btn {
            background: linear-gradient(135deg, #f59e0b, #ef4444);
        }
        .deep-btn:hover {
            background: linear-gradient(135deg, #d97706, #dc2626);
        }
        .search-angle-card {
            background: #fefce8;
            border: 1px solid #fde68a;
            border-radius: 8px;
            padding: 8px 12px;
            font-size: 0.8rem;
        }
    </style>
</head>
<body class="bg-gray-50 min-h-screen text-gray-800">

<nav>
    <div class="nav-container">
        <a href="https://43.99.62.219/" class="back-btn">
            <i class="fa fa-arrow-left"></i> Back to Animal Hormone DB
        </a>
    </div>
</nav>

<main class="max-w-6xl mx-auto px-6 py-8">

    <!-- Search Box -->
    <div class="glass rounded-2xl p-6 mb-8 pulse-glow">
        <div class="flex flex-col gap-4">
            <div class="flex items-center gap-2 text-sm text-gray-500">
                <svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z"/></svg>
                Enter keywords in any language — AI will optimize search terms for PubMed
            </div>
            <div class="flex gap-3">
                <input id="searchInput" type="text"
                    placeholder="e.g. 胰岛素相关基因, insulin resistance, STAT3信号通路..."
                    class="flex-1 bg-white border border-gray-200 rounded-xl px-5 py-3 text-gray-800 placeholder-gray-400 focus:outline-none focus:border-violet-400 focus:ring-2 focus:ring-violet-100 transition-all"
                    onkeydown="if(event.key==='Enter') doSearch('normal')"
                />
                <button onclick="doSearch('normal')" id="searchBtn"
                    class="bg-gradient-to-r from-violet-600 to-blue-600 hover:from-violet-500 hover:to-blue-500 text-white px-6 py-3 rounded-xl font-medium transition-all flex items-center gap-2 whitespace-nowrap">
                    <svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z"/></svg>
                    Search
                </button>
                <button onclick="doSearch('deep')" id="deepSearchBtn"
                    class="deep-btn text-white px-6 py-3 rounded-xl font-medium transition-all flex items-center gap-2 whitespace-nowrap">
                    <svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z"/></svg>
                    Deep Search
                </button>
            </div>
            <div class="flex gap-2 flex-wrap">
                <span class="text-xs text-gray-400">Popular:</span>
                <button onclick="quickSearch('胰岛素相关基因')" class="text-xs px-3 py-1 rounded-full border border-gray-200 bg-white hover:border-violet-300 text-gray-600 transition-all cursor-pointer">胰岛素相关基因</button>
                <button onclick="quickSearch('PCOS biomarkers')" class="text-xs px-3 py-1 rounded-full border border-gray-200 bg-white hover:border-violet-300 text-gray-600 transition-all cursor-pointer">PCOS Biomarkers</button>
                <button onclick="quickSearch('STAT3 signaling pathway')" class="text-xs px-3 py-1 rounded-full border border-gray-200 bg-white hover:border-violet-300 text-gray-600 transition-all cursor-pointer">STAT3 Pathway</button>
                <button onclick="quickSearch('apelin receptor obesity')" class="text-xs px-3 py-1 rounded-full border border-gray-200 bg-white hover:border-violet-300 text-gray-600 transition-all cursor-pointer">Apelin & Obesity</button>
            </div>
        </div>
    </div>

    <!-- Loading -->
    <div id="loading" class="hidden flex flex-col items-center justify-center py-20">
        <div class="loader mb-4"></div>
        <p class="text-gray-500" id="loadingText">Searching PubMed...</p>
        <p class="text-gray-400 text-xs mt-2" id="loadingSubText"></p>
    </div>

    <!-- Results -->
    <div id="results" class="hidden space-y-6">
        
        <!-- Search Info -->
        <div id="searchInfoSection" class="glass rounded-2xl p-6 fade-in">
            <h2 class="text-lg font-semibold text-violet-600 mb-3 flex items-center gap-2">
                <svg class="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z"/></svg>
                Search Strategy
            </h2>
            <div id="searchInfoContent"></div>
        </div>

        <!-- Summary -->
        <div id="summarySection" class="glass rounded-2xl p-6 fade-in">
            <h2 class="text-lg font-semibold text-violet-600 mb-3 flex items-center gap-2">
                <svg class="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z"/></svg>
                AI Summary
            </h2>
            <p id="summaryText" class="text-gray-600 leading-relaxed"></p>
        </div>

        <!-- Key Findings -->
        <div id="findingsSection" class="glass rounded-2xl p-6 fade-in">
            <h2 class="text-lg font-semibold text-blue-600 mb-3 flex items-center gap-2">
                <svg class="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M13 10V3L4 14h7v7l9-11h-7z"/></svg>
                Key Findings
            </h2>
            <ul id="findingsList" class="space-y-2"></ul>
        </div>

        <!-- Gene-Literature Table -->
        <div id="tableSection" class="glass rounded-2xl p-6 fade-in">
            <div class="flex items-center justify-between mb-4">
                <h2 class="text-lg font-semibold text-emerald-600 flex items-center gap-2">
                    <svg class="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M3 10h18M3 14h18m-9-4v8m-7 0h14a2 2 0 002-2V8a2 2 0 00-2-2H5a2 2 0 00-2 2v8a2 2 0 002 2z"/></svg>
                    Gene-Literature Reference Table
                    <span id="geneTableCount" class="text-sm text-gray-400 font-normal"></span>
                </h2>
                <button onclick="downloadCSV()" id="downloadBtn" class="text-xs text-gray-500 border border-gray-200 bg-white hover:border-violet-300 px-3 py-1 rounded-lg transition-all flex items-center gap-1">
                    <svg class="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4"/></svg>
                    Download CSV
                </button>
            </div>
            <div class="overflow-x-auto">
                <table class="w-full text-sm" id="geneTable">
                    <thead>
                        <tr class="border-b-2 border-gray-200">
                            <th class="text-left py-2 px-3 text-violet-600 font-medium">Gene Symbol</th>
                            <th class="text-left py-2 px-3 text-gray-600 font-medium">Full Name</th>
                            <th class="text-left py-2 px-3 text-gray-600 font-medium">Function</th>
                            <th class="text-left py-2 px-3 text-amber-600 font-medium">Supporting Literature</th>
                        </tr>
                    </thead>
                    <tbody id="geneTableBody"></tbody>
                </table>
            </div>
        </div>

        <!-- Genes Detail -->
        <div id="genesSection" class="glass rounded-2xl p-6 fade-in">
            <h2 class="text-lg font-semibold text-emerald-600 mb-4 flex items-center gap-2">
                <svg class="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z"/></svg>
                Gene Details & Evidence Chain
            </h2>
            <div id="genesList" class="space-y-4"></div>
        </div>

        <!-- Articles -->
        <div id="articlesSection" class="glass rounded-2xl p-6 fade-in">
            <div class="flex items-center justify-between mb-4">
                <h2 class="text-lg font-semibold text-amber-600 flex items-center gap-2">
                    <svg class="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M12 6.253v13m0-13C10.832 5.477 9.246 5 7.5 5S4.168 5.477 3 6.253v13C4.168 18.477 5.754 18 7.5 18s3.332.477 4.5 1.253m0-13C13.168 5.477 14.754 5 16.5 5c1.747 0 3.332.477 4.5 1.253v13C19.832 18.477 18.247 18 16.5 18c-1.746 0-3.332.477-4.5 1.253"/></svg>
                    Related Articles
                    <span id="articleCount" class="text-sm text-gray-400 font-normal"></span>
                </h2>
                <button onclick="toggleArticles()" id="toggleBtn" class="text-xs text-gray-500 border border-gray-200 bg-white hover:border-violet-300 px-3 py-1 rounded-lg transition-all">
                    Expand / Collapse
                </button>
            </div>
            <div id="articlesList" class="space-y-3"></div>
        </div>
    </div>

    <!-- Error -->
    <div id="error" class="hidden glass rounded-2xl p-8 text-center">
        <p class="text-red-500 text-lg mb-2">Search Failed</p>
        <p id="errorText" class="text-gray-500"></p>
    </div>

</main>

<footer class="border-t border-gray-200 py-4 mt-12 bg-white">
    <div class="max-w-6xl mx-auto px-6 flex items-center justify-between text-xs text-gray-400">
        <span>GeneLit v2.0 · PubMed + DeepSeek AI</span>
        <span>Inner Mongolia University of Science and Technology</span>
    </div>
</footer>

<script>
let articlesVisible = true;

function quickSearch(query) {
    document.getElementById('searchInput').value = query;
    doSearch('normal');
}

async function doSearch(mode) {
    const query = document.getElementById('searchInput').value.trim();
    if (!query) return;

    const btn = document.getElementById('searchBtn');
    const deepBtn = document.getElementById('deepSearchBtn');
    const loading = document.getElementById('loading');
    const results = document.getElementById('results');
    const error = document.getElementById('error');

    btn.disabled = true; btn.classList.add('opacity-50');
    deepBtn.disabled = true; deepBtn.classList.add('opacity-50');
    loading.classList.remove('hidden');
    results.classList.add('hidden');
    error.classList.add('hidden');

    try {
        if (mode === 'deep') {
            // 深度检索
            document.getElementById('loadingText').textContent = '🧠 AI generating multi-angle search terms...';
            document.getElementById('loadingSubText').textContent = 'Step 1/3: DeepSeek analyzing your query';

            const deepResp = await fetch('/api/deep-search', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({query: query, max_results_per_query: 10})
            });

            document.getElementById('loadingText').textContent = '📚 Searching PubMed with multiple queries...';
            document.getElementById('loadingSubText').textContent = 'Step 2/3: Collecting & deduplicating articles';

            const data = await deepResp.json();
            if (!deepResp.ok) throw new Error(data.detail || 'Deep search failed');

            window._lastAnalysis = data.gene_analysis;
            window._lastData = data;
            window._deepMode = true;

            loading.classList.add('hidden');
            results.classList.remove('hidden');
            renderResults(data);

        } else {
            // 普通检索（已含AI优化检索词）
            document.getElementById('loadingText').textContent = '🧠 AI optimizing search terms...';
            document.getElementById('loadingSubText').textContent = 'Step 1/3: DeepSeek generating optimal PubMed query';

            const searchResp = await fetch(`/api/search?query=${encodeURIComponent(query)}&max_results=15&include_ai=false`);

            document.getElementById('loadingText').textContent = '📚 Searching PubMed...';
            document.getElementById('loadingSubText').textContent = 'Step 2/3: Fetching articles with AI-optimized terms';

            const searchData = await searchResp.json();
            if (!searchResp.ok) throw new Error(searchData.detail || 'Search failed');

            if (!searchData.articles || searchData.articles.length === 0) {
                loading.classList.add('hidden');
                results.classList.remove('hidden');
                window._deepMode = false;
                window._lastData = searchData;
                renderResults(searchData);
                return;
            }

            document.getElementById('loadingText').textContent = '🧬 AI analyzing articles for gene extraction...';
            document.getElementById('loadingSubText').textContent = 'Step 3/3: DeepSeek extracting genes & evidence (may take 1-2 min)';

            const aiResp = await fetch(`/api/search?query=${encodeURIComponent(query)}&max_results=15&include_ai=true`);
            const data = await aiResp.json();

            window._lastAnalysis = data.gene_analysis;
            window._lastData = data;
            window._deepMode = false;

            loading.classList.add('hidden');
            results.classList.remove('hidden');
            renderResults(data);
        }

    } catch (e) {
        loading.classList.add('hidden');
        error.classList.remove('hidden');
        document.getElementById('errorText').textContent = e.message;
    } finally {
        btn.disabled = false; btn.classList.remove('opacity-50');
        deepBtn.disabled = false; deepBtn.classList.remove('opacity-50');
    }
}

function renderResults(data) {
    const analysis = data.gene_analysis;
    const isDeep = window._deepMode;

    // Search Info
    let infoHtml = '';
    if (data.optimized_query) {
        const oq = data.optimized_query;
        infoHtml += `<div class="mb-3"><span class="text-sm text-gray-500">Your query:</span> <span class="font-medium text-gray-700">"${data.query}"</span></div>`;
        infoHtml += `<div class="mb-2"><span class="text-sm text-gray-500">AI translated:</span> <span class="font-medium text-blue-600">"${oq.translated || ''}"</span></div>`;
        infoHtml += `<div class="mb-2"><span class="text-sm text-gray-500">PubMed search term:</span> <code class="text-sm bg-violet-50 text-violet-700 px-2 py-0.5 rounded">${oq.search_term || ''}</code></div>`;
        if (oq.related_terms && oq.related_terms.length > 0) {
            infoHtml += `<div class="mt-2"><span class="text-sm text-gray-500">Related terms:</span><div class="flex flex-wrap gap-1 mt-1">`;
            oq.related_terms.forEach(t => { infoHtml += `<span class="text-xs px-2 py-0.5 rounded bg-blue-50 text-blue-600 border border-blue-200">${t}</span>`; });
            infoHtml += `</div></div>`;
        }
    }
    if (data.search_details) {
        infoHtml += `<div class="mb-2"><span class="text-sm text-gray-500">Your query:</span> <span class="font-medium text-gray-700">"${data.query}"</span></div>`;
        infoHtml += `<div class="mb-2"><span class="text-sm font-medium text-amber-600">🔬 Deep Search — ${data.search_details.length} search angles</span></div>`;
        infoHtml += `<div class="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-2 mt-3">`;
        data.search_details.forEach(sd => {
            infoHtml += `<div class="search-angle-card">
                <div class="font-medium text-amber-700">${sd.label}</div>
                <code class="text-xs text-gray-600">${sd.search_term}</code>
                <div class="text-xs text-gray-400 mt-1">Found: ${sd.found} articles (${sd.new_unique} unique)</div>
            </div>`;
        });
        infoHtml += `</div>`;
    }
    if (infoHtml) {
        document.getElementById('searchInfoContent').innerHTML = infoHtml;
        document.getElementById('searchInfoSection').classList.remove('hidden');
    } else {
        document.getElementById('searchInfoSection').classList.add('hidden');
    }

    // Summary
    if (analysis && analysis.summary) {
        document.getElementById('summaryText').textContent = analysis.summary;
        document.getElementById('summarySection').classList.remove('hidden');
    } else {
        document.getElementById('summarySection').classList.add('hidden');
    }

    // Key Findings
    if (analysis && analysis.key_findings && analysis.key_findings.length > 0) {
        const findingsHtml = analysis.key_findings.map(f =>
            `<li class="flex items-start gap-2 text-gray-600 text-sm"><span class="text-blue-500 mt-1">▸</span><span>${f}</span></li>`
        ).join('');
        document.getElementById('findingsList').innerHTML = findingsHtml;
        document.getElementById('findingsSection').classList.remove('hidden');
    } else {
        document.getElementById('findingsSection').classList.add('hidden');
    }

    // Gene Table
    if (analysis && analysis.genes && analysis.genes.length > 0) {
        const genes = analysis.genes;
        document.getElementById('geneTableCount').textContent = `${genes.length} genes found`;
        const tbodyHtml = genes.map(gene => {
            const evLinks = (gene.evidence || []).map(ev =>
                `<a href="https://pubmed.ncbi.nlm.nih.gov/${ev.pmid}/" target="_blank" class="text-violet-600 hover:text-violet-800 underline decoration-violet-300">PMID:${ev.pmid}</a>`
            ).join('<br>');
            return `<tr class="border-b border-gray-100 hover:bg-violet-50/30 transition-colors">
                <td class="py-2 px-3"><span class="gene-tag px-2 py-0.5 rounded font-mono font-bold text-violet-700 text-xs">${gene.symbol}</span></td>
                <td class="py-2 px-3 text-gray-500 text-xs">${gene.full_name || '-'}</td>
                <td class="py-2 px-3 text-gray-600 text-xs" style="max-width:280px">${gene.role || '-'}</td>
                <td class="py-2 px-3 text-xs leading-relaxed">${evLinks || '-'}</td>
            </tr>`;
        }).join('');
        document.getElementById('geneTableBody').innerHTML = tbodyHtml;
        document.getElementById('tableSection').classList.remove('hidden');
    } else {
        document.getElementById('tableSection').classList.add('hidden');
    }

    // Genes Detail
    if (analysis && analysis.genes && analysis.genes.length > 0) {
        const genesHtml = analysis.genes.map(gene => {
            const evidenceHtml = (gene.evidence || []).map(ev =>
                `<div class="evidence-line pl-3 py-1 ml-2 text-sm">
                    <span class="text-gray-500">${ev.finding}</span>
                    <a href="https://pubmed.ncbi.nlm.nih.gov/${ev.pmid}/" target="_blank" class="text-violet-600 hover:text-violet-800 ml-1 text-xs">PMID:${ev.pmid}</a>
                </div>`
            ).join('');
            return `<div class="glass rounded-xl p-4 glass-hover transition-all">
                <div class="flex items-start justify-between mb-2">
                    <div class="flex items-center gap-3">
                        <span class="gene-tag px-3 py-1 rounded-lg font-mono font-bold text-violet-700">${gene.symbol}</span>
                        <span class="text-gray-500 text-sm">${gene.full_name || ''}</span>
                    </div>
                </div>
                <p class="text-gray-600 text-sm mb-2">${gene.role || ''}</p>
                ${evidenceHtml ? `<div class="mt-2"><span class="text-xs text-gray-400 mb-1 block">Evidence:</span>${evidenceHtml}</div>` : ''}
            </div>`;
        }).join('');
        document.getElementById('genesList').innerHTML = genesHtml;
        document.getElementById('genesSection').classList.remove('hidden');
    } else {
        document.getElementById('genesSection').classList.add('hidden');
    }

    // Articles
    const totalArticles = data.total_unique_articles || data.total_articles || (data.articles ? data.articles.length : 0);
    document.getElementById('articleCount').textContent = `${totalArticles} articles`;
    const articlesHtml = (data.articles || []).map((article, i) =>
        `<div class="glass rounded-xl p-4 glass-hover transition-all">
            <div class="flex items-start justify-between">
                <h3 class="text-sm font-medium text-gray-700 flex-1 mr-3">${i+1}. ${article.title}</h3>
                <a href="${article.url}" target="_blank" class="text-violet-600 hover:text-violet-800 text-xs whitespace-nowrap">PubMed ↗</a>
            </div>
            <div class="flex items-center gap-3 mt-2 text-xs text-gray-400">
                <span>${article.authors}</span><span>·</span>
                <span class="text-amber-600">${article.journal}</span><span>·</span>
                <span>${article.year}</span>
                ${article.doi ? `<span>·</span><span>DOI: ${article.doi}</span>` : ''}
            </div>
            <p class="text-gray-500 text-xs mt-2 line-clamp-2">${article.abstract}</p>
            ${article.mesh_terms && article.mesh_terms.length > 0 ?
                `<div class="flex flex-wrap gap-1 mt-2">
                    ${article.mesh_terms.slice(0,5).map(t => `<span class="text-xs px-2 py-0.5 rounded bg-gray-100 text-gray-500">${t}</span>`).join('')}
                </div>` : ''}
        </div>`
    ).join('');
    document.getElementById('articlesList').innerHTML = articlesHtml;
    document.getElementById('articlesSection').classList.remove('hidden');
}

function downloadCSV() {
    const analysis = window._lastAnalysis;
    if (!analysis || !analysis.genes) return;
    let csv = '\\uFEFF';
    csv += 'Gene Symbol,Full Name,Function,Supporting Literature (PMID),Number of References\\n';
    analysis.genes.forEach(gene => {
        const pmids = (gene.evidence || []).map(ev => ev.pmid).join('; ');
        const refCount = (gene.evidence || []).length;
        const fullName = (gene.full_name || '').replace(/,/g, ';');
        const role = (gene.role || '').replace(/,/g, ';');
        csv += `${gene.symbol},"${fullName}","${role}","${pmids}",${refCount}\\n`;
    });
    const blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    const query = document.getElementById('searchInput').value.trim().replace(/[^a-zA-Z0-9 ]/g, '').replace(/ +/g, '_');
    a.href = url;
    a.download = `GeneLit_${query}_${new Date().toISOString().slice(0,10)}.csv`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
    const btn = document.getElementById('downloadBtn');
    const origHtml = btn.innerHTML;
    btn.innerHTML = '<svg class="w-3 h-3" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M5 13l4 4L19 7"/></svg> Downloaded';
    setTimeout(() => { btn.innerHTML = origHtml; }, 2000);
}

function toggleArticles() {
    const list = document.getElementById('articlesList');
    articlesVisible = !articlesVisible;
    list.style.display = articlesVisible ? 'block' : 'none';
}
</script>
</body>
</html>
"""
if __name__ == "__main__":
    uvicorn.run("app:app", host="0.0.0.0", port=8005, workers=2)
