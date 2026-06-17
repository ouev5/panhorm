#!/usr/bin/env python3
"""
Venn/UpSet Analysis Tool
GO/KEGG enrichment analysis via Enrichr API (academic-grade)
Support gene list file upload
"""

from flask import Flask, render_template, request, jsonify
import json
import requests
from itertools import combinations
import tempfile
import subprocess
import base64
import pymysql
import os
import re
import logging
import sys
sys.path.insert(0, '/www/wwwroot/venn-tool')
from single_cell_multi import query_cells_by_genes, compute_tsne_coords
from urllib.parse import quote

# Logging config
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__)

# Config constants
ENRICHR_API = "https://maayanlab.cloud/Enrichr"
STRING_API = "https://string-db.org/api"
REQUEST_TIMEOUT = 60  # seconds (Enrichr can be slow)
MAX_CONTENT_LENGTH = 10 * 1024 * 1024  # 10 MB
GENE_PATTERN = re.compile(r"^[A-Za-z0-9_\-\.]+$")  # Allow alphanumeric, -, _, .

# App config
app.config['MAX_CONTENT_LENGTH'] = MAX_CONTENT_LENGTH

# Enrichr field index constants

# ========== ç©ç§æ å° ==========
SPECIES_MAP = {
    "human":     {"name": "Human",      "string_id": 9606,  "kegg": "KEGG_2021_Human",                   "go": "GO_Biological_Process_2023",          "emoji": "ð§¬"},
    "mouse":     {"name": "Mouse",      "string_id": 10090, "kegg": "KEGG_2019_Mouse",                   "go": "GO_Biological_Process_2023",          "emoji": "ð­"},
    "rat":       {"name": "Rat",        "string_id": 10116, "kegg": "KEGG_2019_Rat",                     "go": "GO_Biological_Process_2023",          "emoji": "ð"},
    "zebrafish": {"name": "Zebrafish",  "string_id": 7955,  "kegg": "KEGG_2019_Zebrafish_Danio_rerio",  "go": "GO_Biological_Process_2023",          "emoji": "ð"},
    "fly":       {"name": "Fly",        "string_id": 7227,  "kegg": "KEGG_2019_Drosophila_melanogaster", "go": "GO_Biological_Process_2023",          "emoji": "ðª°"},
    "worm":      {"name": "C. elegans", "string_id": 6239,  "kegg": "KEGG_2019_Caenorhabditis_elegans", "go": "GO_Biological_Process_2023",          "emoji": "ðª±"},
}
DEFAULT_SPECIES = "human"


# ========== 动态物种支持 ==========
# 主功能的基础已知物种配置；数据库中的其余物种会动态加入。
# 对于没有手工配置 NCBI Taxon ID / Enrichr KEGG 库的物种，后端会降级使用通用 GO 库，
# STRING/KEGG 等依赖特定物种 ID 的功能会使用该物种名称继续展示，但外部库可能无结果。
_SPECIES_CACHE = {"data": None, "ts": 0}

# 常见数据库物种到内部 key 的别名，避免 Human/Mouse/Rat 等重复出现
ORGANISM_KEY_ALIASES = {
    "Human": "human",
    "Mouse": "mouse",
    "Rat": "rat",
    "Zebrafish": "zebrafish",
    "Fruit fly": "fly",
    "Caenorhabditis elegans": "worm",
}

# 部分常见动物的 NCBI Taxon ID；其余 481 个数据库物种仍会列入下拉框
# 如后续需要提高 STRING 准确性，可继续在这里补充 taxon id。
COMMON_TAXON_IDS = {
    "Human": 9606, "Mouse": 10090, "Rat": 10116, "Zebrafish": 7955, "Fruit fly": 7227,
    "Caenorhabditis elegans": 6239, "Bovine": 9913, "Pig": 9823, "Sheep": 9940,
    "Chicken": 9031, "Dog": 9615, "Rabbit": 9986, "Cat": 9685, "Horse": 9796,
    "Rhesus macaque": 9544, "Crab-eating macaque": 9541, "Guinea pig": 10141,
    "Chimpanzee": 9598, "Goat": 9925, "Golden hamster": 10036, "Common carp": 7962,
    "Rainbow trout": 8022, "Atlantic salmon": 8030, "Nile tilapia": 8128,
    "Japanese rice fish": 8090, "African clawed frog": 8355, "Western clawed frog": 8364,
}

TAXON_CACHE_FILE = "/www/wwwroot/venn-tool/data/species_taxon_ids.json"
_TAXON_FILE_CACHE = {"mtime": 0, "data": {}}


def load_taxon_id_cache():
    """读取联网补全的 NCBI Taxonomy ID 缓存文件。"""
    try:
        if not os.path.exists(TAXON_CACHE_FILE):
            return {}
        mtime = os.path.getmtime(TAXON_CACHE_FILE)
        if _TAXON_FILE_CACHE["data"] and _TAXON_FILE_CACHE["mtime"] == mtime:
            return _TAXON_FILE_CACHE["data"]
        with open(TAXON_CACHE_FILE, "r", encoding="utf-8") as f:
            raw = json.load(f)
        data = {}
        for name, info in raw.items():
            taxid = info.get("taxid") if isinstance(info, dict) else info
            if taxid:
                try:
                    data[name] = int(taxid)
                except Exception:
                    pass
        _TAXON_FILE_CACHE["mtime"] = mtime
        _TAXON_FILE_CACHE["data"] = data
        return data
    except Exception as e:
        logger.warning(f"Load taxon id cache failed: {e}")
        return {}


def _slugify_species(name):
    """把数据库 organism 名称转换成前端 option value。"""
    key = re.sub(r"[^A-Za-z0-9]+", "_", (name or "").strip().lower()).strip("_")
    return key or "unknown"


def get_dynamic_species_map(force=False):
    """返回主功能支持的物种映射：基础物种 + hormone_data 数据库全部 organism。"""
    import time
    now = time.time()
    if (not force) and _SPECIES_CACHE["data"] is not None and now - _SPECIES_CACHE["ts"] < 300:
        return _SPECIES_CACHE["data"]

    dynamic = {k: dict(v) for k, v in SPECIES_MAP.items()}
    taxon_cache = load_taxon_id_cache()
    try:
        conn = get_hormone_db_connection()
        if conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT organism, COUNT(*) as cnt
                FROM hormone_data
                WHERE organism IS NOT NULL AND organism != ''
                GROUP BY organism
                ORDER BY cnt DESC
            """)
            for row in cursor.fetchall():
                org = row.get("organism")
                if not org:
                    continue
                key = ORGANISM_KEY_ALIASES.get(org) or _slugify_species(org)
                # 避免 key 冲突
                base_key = key
                n = 2
                while key in dynamic and dynamic[key].get("db_name") not in (None, org) and dynamic[key].get("name") != org:
                    key = f"{base_key}_{n}"
                    n += 1
                if key in dynamic:
                    dynamic[key]["db_count"] = row.get("cnt", 0)
                    dynamic[key]["db_name"] = org
                    dynamic[key]["from_db"] = True
                    if not dynamic[key].get("string_id") and taxon_cache.get(org):
                        dynamic[key]["string_id"] = taxon_cache.get(org)
                    continue
                taxid = COMMON_TAXON_IDS.get(org) or taxon_cache.get(org)
                dynamic[key] = {
                    "name": org,
                    "db_name": org,
                    "db_count": row.get("cnt", 0),
                    "from_db": True,
                    "string_id": taxid,  # None 表示暂未配置 Taxon ID
                    "kegg": "KEGG_2021_Human",  # Enrichr 不提供 481 个物种专属库，未知物种用兜底库避免接口崩溃
                    "go": "GO_Biological_Process_2023",
                    "emoji": "🧬",
                }
            cursor.close()
            conn.close()
    except Exception as e:
        logger.warning(f"Load dynamic species from DB failed: {e}")

    _SPECIES_CACHE["data"] = dynamic
    _SPECIES_CACHE["ts"] = now
    return dynamic


def get_species_info(species_key):
    """获取物种配置，优先动态数据库物种，找不到时回到默认 human。"""
    return get_dynamic_species_map().get(species_key, get_dynamic_species_map().get(DEFAULT_SPECIES, SPECIES_MAP[DEFAULT_SPECIES]))

ENRICHR_IDX = {
    "index": 0,
    "term": 1,
    "pvalue": 2,
    "z_score": 3,
    "combined_score": 4,
    "genes": 5,
    "adj_pvalue": 6,
    "old_pvalue": 7,
    "old_adj_pvalue": 8
}


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/parse_file", methods=["POST"])
def parse_file():
    """Parse uploaded gene file"""
    try:
        if "file" not in request.files:
            return jsonify({"error": "No file uploaded"}), 400

        file = request.files["file"]
        if file.filename == "":
            return jsonify({"error": "No file selected"}), 400

        # Read file content
        content = file.read().decode("utf-8", errors="ignore")

        # Parse gene names
        genes = []
        lines = content.strip().split("\n")

        for line in lines:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # Handle CSV format (comma separated)
            if "," in line:
                parts = line.split(",")
                for part in parts:
                    gene = part.strip().strip("\"'")
                    if gene and not gene.startswith("#"):
                        genes.append(gene)
            # Handle tab separated
            elif "\t" in line:
                parts = line.split("\t")
                gene = parts[0].strip().strip("\"'")
                if gene and not gene.startswith("#"):
                    genes.append(gene)
            # Handle one gene per line
            else:
                gene = line.strip().strip("\"'")
                if gene and not gene.startswith("#"):
                    genes.append(gene)

        # Deduplicate and filter invalid gene names
        genes = list(dict.fromkeys(genes))
        genes = [g for g in genes if GENE_PATTERN.match(g)]

        logger.info(f"File parsed, valid genes: {len(genes)}")
        return jsonify({"success": True, "genes": genes, "count": len(genes)})

    except Exception as e:
        logger.error(f"File parsing failed: {str(e)}")
        return jsonify({"error": f"File parsing failed: {str(e)}"}), 500


@app.route("/api/analyze", methods=["POST"])
def analyze():
    """Analyze multiple lists, compute intersections and unions"""
    try:
        data = request.get_json()
        lists = data.get("lists", {})

        sets = {name: set(items) for name, items in lists.items() if items}

        if len(sets) < 2:
            return jsonify({"error": "At least 2 lists required"}), 400

        result = {
            "lists": {name: list(items) for name, items in sets.items()},
            "sizes": {name: len(items) for name, items in sets.items()},
            # Exclusive UpSet intersections: genes belonging to exactly this set combination.
            "intersections": [],
            # Inclusive intersections: genes shared by this combination, regardless of other lists.
            # This is useful when 4+ lists are highly overlapping; otherwise the classic exclusive
            # UpSet plot can collapse to a single all-lists column.
            "intersections_inclusive": [],
            "union": list(set.union(*sets.values())) if sets else []
        }

        list_names = list(sets.keys())

        for r in range(1, len(list_names) + 1):
            for combo in combinations(list_names, r):
                intersection = set.intersection(*[sets[n] for n in combo])
                if intersection:
                    result["intersections_inclusive"].append({
                        "sets": list(combo),
                        "size": len(intersection),
                        "elements": sorted(intersection)
                    })

                others = set.union(*[sets[n] for n in list_names if n not in combo])                     if len(combo) < len(list_names) else set()
                exclusive = intersection - others

                if exclusive:
                    result["intersections"].append({
                        "sets": list(combo),
                        "size": len(exclusive),
                        "elements": sorted(exclusive)
                    })

        result["intersections"].sort(key=lambda x: x["size"], reverse=True)
        result["intersections_inclusive"].sort(key=lambda x: (x["size"], len(x["sets"])), reverse=True)
        logger.info(
            f"Intersection analysis done, exclusive: {len(result['intersections'])}, "
            f"inclusive: {len(result['intersections_inclusive'])}"
        )
        return jsonify(result)

    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        return jsonify({"error": f"Analysis failed: {str(e)}"}), 500




@app.route("/api/species", methods=["GET"])
def get_species():
    """Return supported species list: built-in species + all organisms from hormone_data DB"""
    result = []
    dyn = get_dynamic_species_map(force=True)
    for key, info in dyn.items():
        result.append({
            "key": key,
            "name": info.get("name", key),
            "db_name": info.get("db_name", info.get("name", key)),
            "emoji": info.get("emoji", "🧬"),
            "string_id": info.get("string_id"),
            "count": info.get("db_count", 0),
            "from_db": info.get("from_db", False),
        })
    # 数据库物种优先按 count 降序，基础物种没有 count 时排在前面
    result.sort(key=lambda x: (0 if x.get("count", 0) == 0 else 1, -int(x.get("count", 0) or 0), x.get("name", "")))
    return jsonify({"species": result, "default": DEFAULT_SPECIES, "total": len(result)})

@app.route("/api/enrichment", methods=["POST"])
def enrichment():
    """GO/KEGG enrichment via Enrichr"""
    try:
        data = request.get_json()
        genes = data.get("genes", [])
        analysis_type = data.get("type", "go")
        species_key = data.get("species", DEFAULT_SPECIES)

        if not genes:
            return jsonify({"error": "Please provide a gene list"}), 400

        if len(genes) > 1000:
            logger.warning(f"Gene count {len(genes)} >1000, Enrichr may truncate")

        # Select gene library by species
        sp = get_species_info(species_key)
        if analysis_type == "go":
            gene_library = sp["go"]
        elif analysis_type == "kegg":
            gene_library = sp["kegg"]
        else:
            return jsonify({"error": "Unsupported enrichment type, use 'go' or 'kegg'"}), 400

        genes_str = "\n".join(genes)
        # Submit gene list with retries
        max_retries = 3
        response = None
        for attempt in range(max_retries):
            try:
                response = requests.post(
                    f"{ENRICHR_API}/addList",
                    files={"list": (None, genes_str), "description": (None, "gene list")},
                    timeout=REQUEST_TIMEOUT
                )
                if response.status_code == 200:
                    break
                logger.warning(f"Enrichr submit attempt {attempt+1} failed, status: {response.status_code}")
            except requests.exceptions.Timeout:
                logger.warning(f"Enrichr submit attempt {attempt+1} timed out")
            if attempt < max_retries - 1:
                import time; time.sleep(3)

        if not response or response.status_code != 200:
            logger.error(f"Enrichr submit failed after {max_retries} retries, status: {response.status_code if response else 'N/A'}")
            return jsonify({"error": f"Failed to submit gene list (retried {max_retries} times), Enrichr service unavailable, please try later"}), 502

        user_list_id = response.json().get("userListId")
        if not user_list_id:
            return jsonify({"error": "Failed to get user list ID"}), 500

        response = requests.get(
            f"{ENRICHR_API}/enrich",
            params={"userListId": user_list_id, "backgroundType": gene_library},
            timeout=REQUEST_TIMEOUT
        )

        if response.status_code != 200:
            logger.error(f"Enrichr get results failed, status: {response.status_code}")
            return jsonify({"error": f"Failed to get enrichment results, status: {response.status_code}"}), 500

        enrichr_result = response.json()
        results = enrichr_result.get(gene_library, [])

        enriched = []
        for item in results:
            # Safely get fields
            term = item[ENRICHR_IDX["term"]] if len(item) > ENRICHR_IDX["term"] else "N/A"
            pvalue = item[ENRICHR_IDX["pvalue"]] if len(item) > ENRICHR_IDX["pvalue"] else 1.0
            adj_pvalue = item[ENRICHR_IDX["adj_pvalue"]] if len(item) > ENRICHR_IDX["adj_pvalue"] else pvalue
            zscore = round(item[ENRICHR_IDX["z_score"]], 2) if len(item) > ENRICHR_IDX["z_score"] else 0.0
            combined_score = round(item[ENRICHR_IDX["combined_score"]], 1) if len(item) > ENRICHR_IDX["combined_score"] else 0.0
            genes_overlap = item[ENRICHR_IDX["genes"]] if len(item) > ENRICHR_IDX["genes"] else []

            enriched.append({
                "term": term,
                "pvalue": pvalue,
                "adj_pvalue": adj_pvalue,
                "zscore": zscore,
                "combined_score": combined_score,
                "intersection_size": len(genes_overlap),
                "genes": genes_overlap
            })

        logger.info(f"Enrichment done, returning {len(enriched)} results")
        return jsonify({"success": True, "results": enriched, "source": "Enrichr"})

    except requests.exceptions.Timeout:
        logger.error("Enrichr request timed out")
        return jsonify({"error": "Enrichr request timed out, please try again"}), 504
    except requests.exceptions.RequestException as e:
        logger.error(f"Enrichr network error: {str(e)}")
        return jsonify({"error": f"Network error: {str(e)}"}), 500
    except Exception as e:
        logger.error(f"Enrichment analysis failed: {str(e)}")
        return jsonify({"error": f"Enrichment analysis failed: {str(e)}"}), 500


@app.route("/api/string_url", methods=["POST"])
def string_url():
    """Generate STRING database URL"""
    try:
        data = request.get_json()
        genes = data.get("genes", [])

        if not genes:
            return jsonify({"error": "Please provide a gene list"}), 400

        # STRING API recommends identifiers param with encoding
        species_key = data.get("species", DEFAULT_SPECIES)
        sp = get_species_info(species_key)
        gene_str = "%0d".join(genes[:50])  # STRING multi-gene separated by %0d (newline) for stability
        string_url = f"https://string-db.org/cgi/network?taskId=&identifiers={gene_str}&species={sp.get('string_id') or 9606}"

        logger.info(f"Generated STRING URL, genes: {len(genes[:50])}, species: {sp['name']}")
        return jsonify({"url": string_url, "genes": genes[:50]})

    except Exception as e:
        logger.error(f"Failed to generate STRING URL: {str(e)}")
        return jsonify({"error": str(e)}), 500


# ========== STRING Network Data API ==========
@app.route("/api/string_network", methods=["POST"])
def string_network():
    """Get protein interaction network from STRING"""
    try:
        data = request.get_json()
        genes = data.get("genes", [])
        min_score = data.get("min_score", 400)  # Minimum confidence score (0-1000)
        species_key = data.get("species", DEFAULT_SPECIES)
        sp = get_species_info(species_key)
        
        if not genes:
            return jsonify({"error": "Please provide a gene list"}), 400
        
        # Limit gene count
        genes = genes[:50]
        
        # Call STRING API for network data
        # 1. Map genes to STRING IDs
        string_ids = []
        for gene in genes:
            try:
                # Get STRING ID
                url = f"{STRING_API}/json/resolve?identifier={quote(gene)}&species={sp.get('string_id') or 9606}"
                resp = requests.get(url, timeout=10)
                if resp.status_code == 200 and resp.json():
                    result = resp.json()[0]
                    string_ids.append({
                        "query": gene,
                        "stringId": result.get("stringId", ""),
                        "preferredName": result.get("preferredName", gene),
                        "ncbiTaxonId": result.get("ncbiTaxonId", 9606)
                    })
            except Exception as e:
                logger.warning(f"STRING ID mapping failed for {gene}: {e}")
                continue
        
        if not string_ids:
            return jsonify({"error": "Unable to map any genes to STRING database"}), 400
        
        # 2. Get interaction data(using interactions endpoint for within-list gene interactions)
        interactions = []
        preferred_names = [s["preferredName"] for s in string_ids]
        
        # Use interactions endpoint: query interactions within gene list
        try:
            gene_str = "%0d".join(preferred_names)
            url = f"{STRING_API}/json/interactions?identifiers={gene_str}&species={sp.get('string_id') or 9606}&required_score={min_score}"
            resp = requests.get(url, timeout=60)
            if resp.status_code == 200:
                for item in resp.json():
                    interactions.append({
                        "stringId_A": item.get("stringId_A", ""),
                        "stringId_B": item.get("stringId_B", ""),
                        "preferredName_A": item.get("preferredName_A", ""),
                        "preferredName_B": item.get("preferredName_B", ""),
                        "score": item.get("score", 0),
                    })
                logger.info(f"STRING interactions: {len(interactions)} records")
        except Exception as e:
            logger.warning(f"STRING interaction fetch failed: {e}")
        
        # 3. Build network data structure (for Cytoscape.js)
        nodes = []
        edges = []
        node_set = set()
        
        # Add nodes
        for s in string_ids:
            if s["stringId"] not in node_set:
                nodes.append({
                    "data": {
                        "id": s["stringId"],
                        "label": s["preferredName"],
                        "query": s["query"]
                    }
                })
                node_set.add(s["stringId"])
        
        # Add edges (deduplicate)
        edge_set = set()
        for inter in interactions:
            id_a = inter["stringId_A"]
            id_b = inter["stringId_B"]
            
            # Ensure both nodes are in our list
            if id_a not in node_set or id_b not in node_set:
                continue
            
            # Create unique edge ID
            edge_key = tuple(sorted([id_a, id_b]))
            if edge_key in edge_set:
                continue
            edge_set.add(edge_key)
            
            edges.append({
                "data": {
                    "id": f"{id_a}_{id_b}",
                    "source": id_a,
                    "target": id_b,
                    "weight": inter["score"] / 1000,  # Normalize to 0-1
                    "score": inter["score"]
                }
            })
        
        logger.info(f"STRING network: {len(nodes)} nodes, {len(edges)} edges")
        
        return jsonify({
            "success": True,
            "nodes": nodes,
            "edges": edges,
            "string_ids": string_ids,
            "stats": {
                "n_genes": len(genes),
                "n_mapped": len(string_ids),
                "n_nodes": len(nodes),
                "n_edges": len(edges)
            }
        })
        
    except requests.exceptions.Timeout:
        logger.error("STRING API request timed out")
        return jsonify({"error": "STRING API request timed out, please try again"}), 504
    except Exception as e:
        logger.error(f"Failed to get STRING network data: {str(e)}")
        return jsonify({"error": str(e)}), 500


@app.route("/api/string_network_tsv", methods=["POST"])
def string_network_tsv():
    """Get TSV format network data from STRING (for import)"""
    try:
        data = request.get_json()
        genes = data.get("genes", [])
        min_score = data.get("min_score", 400)
        
        if not genes:
            return jsonify({"error": "Please provide a gene list"}), 400
        
        genes = genes[:50]
        gene_str = "|".join(quote(g) for g in genes)
        
        # Get TSV format interaction data
        species_key = data.get("species", DEFAULT_SPECIES)
        sp = get_species_info(species_key)
        url = f"{STRING_API}/tsv/interactions?identifiers={gene_str}&species={sp.get('string_id') or 9606}&required_score={min_score}"
        resp = requests.get(url, timeout=30)
        
        if resp.status_code != 200:
            return jsonify({"error": f"STRING API error: {resp.status_code}"}), 500
        
        # Parse TSV
        lines = resp.text.strip().split("\n")
        if len(lines) < 2:
            return jsonify({"success": True, "nodes": [], "edges": [], "raw_tsv": resp.text})
        
        header = lines[0].split("\t")
        nodes = {}
        edges = []
        
        for line in lines[1:]:
            cols = line.split("\t")
            if len(cols) < 10:
                continue
            
            # Extract node and edge info
            name_a = cols[0] if len(cols) > 0 else ""
            name_b = cols[1] if len(cols) > 1 else ""
            score = float(cols[9]) if len(cols) > 9 else 0
            
            # Add nodes
            if name_a and name_a not in nodes:
                nodes[name_a] = {"data": {"id": name_a, "label": name_a}}
            if name_b and name_b not in nodes:
                nodes[name_b] = {"data": {"id": name_b, "label": name_b}}
            
            # Add edges
            if name_a and name_b:
                edges.append({
                    "data": {
                        "id": f"{name_a}_{name_b}",
                        "source": name_a,
                        "target": name_b,
                        "weight": score / 1000,
                        "score": score
                    }
                })
        
        return jsonify({
            "success": True,
            "nodes": list(nodes.values()),
            "edges": edges,
            "stats": {
                "n_nodes": len(nodes),
                "n_edges": len(edges)
            }
        })
        
    except Exception as e:
        logger.error(f"Failed to get STRING TSV data: {str(e)}")
        return jsonify({"error": str(e)}), 500


# ========== Database Connection ==========
def get_hormone_db_connection():
    """Connect to animal hormone database"""
    try:
        conn = pymysql.connect(
            host='127.0.0.1',
            user='hhh',
            password='YOUR_VALUE_HERE',
            database='animal_hormone',
            charset='utf8mb4',
            cursorclass=pymysql.cursors.DictCursor
        )
        return conn
    except Exception as e:
        logger.error(f"Database connection failed: {str(e)}")
        return None

# ========== R Sankey Bubble Chart Generation ==========
@app.route("/api/sankey_bubble_plot", methods=["POST"])
def sankey_bubble_plot():
    """Generate sankey bubble chart using R"""
    try:
        data = request.get_json()
        results = data.get("results", [])
        sort_by = data.get("sort_by", "pvalue")
        sort_order = data.get("sort_order", "asc")
        top_n = data.get("top_n", 25)
        try:
            top_n = max(1, min(100, int(top_n)))
        except Exception:
            top_n = 25
        
        if not results:
            return jsonify({"error": "No enrichment result data"}), 400
        
        r_data = {
            "terms": [r.get("term", "") for r in results],
            "genes": [r.get("genes", []) for r in results],
            "pvalue": [r.get("pvalue", 1) for r in results],
            "adj_pvalue": [r.get("adj_pvalue", r.get("pvalue", 1)) for r in results],
            "combined_score": [r.get("combined_score", 0) for r in results],
            "intersection_size": [r.get("intersection_size", len(r.get("genes", []))) for r in results],
            "sort_by": sort_by,
            "sort_order": sort_order,
            "top_n": top_n,
            "ratio": []
        }
        
        total_genes = data.get("total_genes", 100)
        for r in results:
            gene_count = r.get("intersection_size", len(r.get("genes", [])))
            r_data["ratio"].append(gene_count / max(1, total_genes))
        
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(r_data, f)
            input_file = f.name
        
        output_file = tempfile.mktemp(suffix=".png")
        
        try:
            r_script = "/www/wwwroot/venn-tool/r_scripts/sankey_bubble.R"
            result = subprocess.run(
                ["/root/miniconda3/bin/Rscript", r_script, input_file, output_file],
                capture_output=True, text=True, timeout=180
            )
            
            if result.returncode != 0:
                return jsonify({"error": f"R script execution failed: {result.stderr[:200]}"}), 500
            
            if os.path.exists(output_file):
                with open(output_file, "rb") as f:
                    img_data = base64.b64encode(f.read()).decode("utf-8")
                os.unlink(input_file)
                os.unlink(output_file)
                return jsonify({"success": True, "image": f"data:image/png;base64,{img_data}"})
            else:
                return jsonify({"error": "Image generation failed"}), 500
                
        except subprocess.TimeoutExpired:
            return jsonify({"error": "R script execution timed out"}), 504
        except Exception as e:
            return jsonify({"error": str(e)}), 500
            
    except Exception as e:
        return jsonify({"error": str(e)}), 500


# ========== Database Search API ==========
@app.route("/api/species_list", methods=["GET"])
def species_list():
    """Get all available species from database"""
    try:
        conn = get_hormone_db_connection()
        if not conn:
            return jsonify({"error": "Database connection failed", "species": []}), 500
        cursor = conn.cursor()
        cursor.execute("""
            SELECT DISTINCT organism, COUNT(*) as cnt 
            FROM hormone_data 
            WHERE organism IS NOT NULL AND organism != '' 
            GROUP BY organism 
            ORDER BY cnt DESC
        """)
        species = []
        for row in cursor.fetchall():
            species.append({"name": row["organism"], "count": row["cnt"]})
        cursor.close()
        conn.close()
        return jsonify({"species": species})
    except Exception as e:
        logger.error(f"Species list failed: {str(e)}")
        return jsonify({"error": str(e), "species": []}), 500


@app.route("/api/search_database", methods=["GET"])
def search_database():
    """Search hormones and genes in database"""
    try:
        query = request.args.get("q", "").strip()
        search_type = request.args.get("type", "all")
        organism = request.args.get("organism", "").strip()
        
        if not query:
            return jsonify({"error": "Please enter a search keyword", "results": []}), 400
        
        conn = get_hormone_db_connection()
        if not conn:
            return jsonify({"error": "Database connection failed", "results": []}), 500
        
        cursor = conn.cursor()
        results = []
        
        # Search hormones - deduplicate by (hormone_name, organism)
        if search_type in ["all", "hormone"]:
            if organism:
                cursor.execute("""
                    SELECT hormone_name, organism,
                           GROUP_CONCAT(DISTINCT related_genes) as all_genes
                    FROM hormone_data 
                    WHERE hormone_name LIKE %s AND organism = %s
                    GROUP BY hormone_name, organism
                    ORDER BY hormone_name
                    LIMIT 200
                """, (f"%{query}%", organism))
            else:
                cursor.execute("""
                    SELECT hormone_name, organism,
                           GROUP_CONCAT(DISTINCT related_genes) as all_genes
                    FROM hormone_data 
                    WHERE hormone_name LIKE %s
                    GROUP BY hormone_name, organism
                    ORDER BY hormone_name
                    LIMIT 200
                """, (f"%{query}%",))
            
            for row in cursor.fetchall():
                gene_set = set()
                for g in str(row["all_genes"]).split(","):
                    g = g.strip()
                    if g:
                        gene_set.add(g)
                genes = sorted(gene_set)
                org = row.get("organism", "")
                org_label = f" [{org}]" if org else ""
                results.append({
                    "type": "hormone",
                    "name": row["hormone_name"],
                    "genes": genes[:5],
                    "all_genes": genes,
                    "gene_count": len(genes),
                    "display": f"Hormone: {row['hormone_name']} ({len(genes)} genes){org_label}",
                    "organism": org
                })
        
        # Search genes - deduplicate by (hormone_name, organism, gene)
        if search_type in ["all", "gene"]:
            if organism:
                cursor.execute("""
                    SELECT hormone_name, related_genes, organism
                    FROM hormone_data 
                    WHERE related_genes LIKE %s AND organism = %s
                    GROUP BY hormone_name, organism
                    LIMIT 200
                """, (f"%{query}%", organism))
            else:
                cursor.execute("""
                    SELECT hormone_name, related_genes, organism
                    FROM hormone_data 
                    WHERE related_genes LIKE %s
                    GROUP BY hormone_name, organism
                    LIMIT 200
                """, (f"%{query}%",))
            
            for row in cursor.fetchall():
                genes = [g.strip() for g in str(row["related_genes"]).split(",") if g.strip()]
                if query.upper() in [g.upper() for g in genes]:
                    org = row.get("organism", "")
                    org_label = f" [{org}]" if org else ""
                    results.append({
                        "type": "gene",
                        "name": query,
                        "hormone": row["hormone_name"],
                        "genes": [query],
                        "display": f"Gene: {query} (from {row['hormone_name']}){org_label}",
                        "organism": org
                    })
        
        cursor.close()
        conn.close()
        
        return jsonify({
            "success": True,
            "query": query,
            "count": len(results),
            "results": results  
        })
        
    except Exception as e:
        logger.error(f"Database search failed: {str(e)}")
        return jsonify({"error": str(e), "results": []}), 500

PBMC_DATA = None

def load_pbmc_data():
    """Load PBMC3k data"""
    global PBMC_DATA
    if PBMC_DATA is not None:
        return PBMC_DATA
    
    import os
    import json
    
    data_file = "/www/wwwroot/venn-tool/data/pbmc3k_umap.json"
    if os.path.exists(data_file):
        with open(data_file, "r") as f:
            PBMC_DATA = json.load(f)
        return PBMC_DATA
    
    # If no data files, generate simulated data
    import random
    import math
    
    cell_types = ["CD4 T cells", "CD8 T cells", "B cells", "NK cells", "Monocytes", "DC", "Platelets"]
    genes_in_data = ["CD3D", "CD3E", "CD4", "CD8A", "CD8B", "CD19", "CD79A", "NKG7", "GNLY", 
                     "CD14", "LYZ", "FCGR3A", "CD1C", "PPBP", "MS4A1", "IL7R", "CCR7",
                     "S100A4", "LEF1", "NKG7", "GZMB", "CST7", "FCER1A", "LGALS2"]
    
    cells = []
    for i in range(2700):
        cell_type = random.choice(cell_types)
        angle = random.uniform(0, 2 * math.pi)
        radius = random.uniform(2, 8)
        
        # Adjust positions by cell type
        type_offset = cell_types.index(cell_type) * (2 * math.pi / len(cell_types))
        x = radius * math.cos(angle + type_offset) + random.gauss(0, 0.5)
        y = radius * math.sin(angle + type_offset) + random.gauss(0, 0.5)
        
        cell = {
            "x": round(x, 4),
            "y": round(y, 4),
            "cell_type": cell_type,
            "genes": [g for g in random.sample(genes_in_data, random.randint(3, 8))]
        }
        cells.append(cell)
    
    PBMC_DATA = {
        "cells": cells,
        "genes": genes_in_data,
        "cell_types": cell_types
    }
    return PBMC_DATA

# ========== Single Cell Analysis API ==========
# Tabula Sapiens dataset (human whole-body single cell atlas)
TS_DATA = None
TS_GENE_MAPPING = None

def load_tabula_sapiens():
    """Load Tabula Sapiens data"""
    global TS_DATA, TS_GENE_MAPPING
    if TS_DATA is not None:
        return TS_DATA, TS_GENE_MAPPING
    
    import os
    import json
    
    data_dir = "/www/wwwroot/venn-tool/single_cell_data"
    
    # Load cell metadata
    meta_file = os.path.join(data_dir, "cell_metadata.json")
    gene_file = os.path.join(data_dir, "gene_list.json")
    
    if os.path.exists(meta_file) and os.path.exists(gene_file):
        with open(meta_file, "r") as f:
            meta = json.load(f)
        with open(gene_file, "r") as f:
            genes = json.load(f)
        # Cell data in meta["cells"]
        TS_DATA = meta.get("cells", [])
        # Create gene index mapping
        TS_GENE_MAPPING = {g.upper(): i for i, g in enumerate(genes)}
        logger.info(f"Loaded Tabula Sapiens: {len(TS_DATA)} cells, {len(genes)} genes")
        return TS_DATA, TS_GENE_MAPPING
    
    # If no data, return empty
    logger.warning("Tabula Sapiens data not found")
    return [], {}



@app.route("/api/single_cell_background", methods=["POST"])
def single_cell_background():
    """Return background UMAP coordinates for all cells (no gene filtering)"""
    try:
        data = request.get_json()
        species_key = data.get("species", DEFAULT_SPECIES)
        
        from single_cell_multi import load_data
        matrix, gene_index, cell_meta = load_data(species_key)
        
        if matrix is None:
            return jsonify({"error": f"No single cell data for species: {species_key}"}), 404
        
        # Return all cells with just coordinates + metadata (no expression)
        # Sample if too many (keep under 100k for performance)
        max_bg = 100000
        n_cells = len(cell_meta)
        
        import random
        if n_cells > max_bg:
            indices = sorted(random.sample(range(n_cells), max_bg))
        else:
            indices = list(range(n_cells))
        
        bg_cells = []
        for i in indices:
            c = cell_meta[i]
            bg_cells.append({
                "cell_id": c.get("cell_id", ""),
                "cell_type": c.get("cell_type", ""),
                "tissue": c.get("tissue", ""),
                "x": c.get("x", 0),
                "y": c.get("y", 0),
            })
        
        all_cell_types = sorted(set(c["cell_type"] for c in bg_cells))
        all_tissues = sorted(set(c["tissue"] for c in bg_cells))
        
        return jsonify({
            "cells": bg_cells,
            "n_cells": len(bg_cells),
            "total_cells": n_cells,
            "cell_types": all_cell_types,
            "tissues": all_tissues,
            "species": species_key
        })
    
    except Exception as e:
        logger.error(f"Background query failed: {str(e)}")
        return jsonify({"error": str(e)}), 500

@app.route("/api/single_cell_query", methods=["POST"])
def single_cell_query():
    """Single cell UMAP query - using real Tabula Sapiens expression data"""
    try:
        data = request.get_json()
        query_genes = data.get("genes", [])
        species_key = data.get("species", DEFAULT_SPECIES)
        
        if not query_genes:
            return jsonify({"error": "Please provide a gene list"}), 400
        
        # Query with real expression matrix (multi-species support)
        cells, found_genes, not_found_genes = query_cells_by_genes(query_genes, species=species_key, max_cells=5000)
        
        if not cells:
            return jsonify({
                "cells": [],
                "n_cells": 0,
                "found_genes": found_genes,
                "not_found_genes": not_found_genes,
                "cell_types": []
            })
        
        # Get all cell types
        all_cell_types = list(set(c.get("cell_type", "unknown") for c in cells))
        
        return jsonify({
            "cells": cells,
            "n_cells": len(cells),
            "found_genes": found_genes,
            "not_found_genes": not_found_genes,
            "cell_types": all_cell_types
        })
        
    except Exception as e:
        logger.error(f"Single cell query failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500

@app.route("/api/single_cell_query_tsne", methods=["POST"])
def single_cell_query_tsne():
    """Single cell t-SNE query - real computation"""
    try:
        data = request.get_json()
        query_genes = data.get("genes", [])
        species_key = data.get("species", DEFAULT_SPECIES)
        
        if not query_genes:
            return jsonify({"error": "Please provide a gene list"}), 400
        
        # Query with real expression matrix
        cells, found_genes, not_found_genes = query_cells_by_genes(query_genes, max_cells=3000)
        
        if not cells:
            return jsonify({
                "cells": [],
                "n_cells": 0,
                "found_genes": found_genes,
                "not_found_genes": not_found_genes,
                "cell_types": []
            })
        
        # Compute real t-SNE coordinates
        cells = compute_tsne_coords(cells)
        
        # Get all cell types
        all_cell_types = list(set(c.get("cell_type", "unknown") for c in cells))
        
        return jsonify({
            "cells": cells,
            "n_cells": len(cells),
            "found_genes": found_genes,
            "not_found_genes": not_found_genes,
            "cell_types": all_cell_types
        })
        
    except Exception as e:
        logger.error(f"t-SNE query failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500



# ========== AI Summary (DeepSeek) ==========
DEEPSEEK_API_KEY = "YOUR_VALUE_HERE"
DEEPSEEK_API_URL = "https://api.deepseek.com/v1/chat/completions"
DEEPSEEK_MODEL = "deepseek-v4-pro"


@app.route("/api/ai_summary", methods=["POST"])
def ai_summary():
    try:
        data = request.get_json(force=True)
        analysis_data = data.get("analysis_data", {})
        enrichment_data = data.get("enrichment_data", {})
        selected_genes = data.get("selected_genes", [])
        species = data.get("species", "human")

        context_parts = []
        context_parts.append(f"Species: {species}")

        if analysis_data:
            sizes = analysis_data.get("sizes", {})
            if sizes:
                for name, count in sizes.items():
                    context_parts.append(f'Gene List "{name}": {count} genes')
            intersections = analysis_data.get("intersections", [])
            if intersections:
                context_parts.append(f"Intersection Analysis ({len(intersections)} intersections):")
                for i, inter in enumerate(intersections[:10]):
                    sets_str = " ∩ ".join(inter.get("sets", []))
                    context_parts.append(f"  {sets_str}: {inter.get('size', 0)} shared genes")

        if selected_genes:
            gene_str = ", ".join(selected_genes[:50])
            if len(selected_genes) > 50:
                gene_str += f" ... (+{len(selected_genes) - 50} more)"
            context_parts.append(f"Currently selected genes ({len(selected_genes)}): {gene_str}")

        if enrichment_data:
            go_results = enrichment_data.get("go", [])
            kegg_results = enrichment_data.get("kegg", [])

            if go_results:
                context_parts.append(f"GO Enrichment ({len(go_results)} terms):")
                for r in go_results[:8]:
                    genes_in = r.get("genes", [])
                    gene_list_str = ", ".join(genes_in[:5]) if isinstance(genes_in, list) else str(genes_in)[:60]
                    context_parts.append(
                        f"  - {r.get('term', 'N/A')} (p={r.get('pvalue', 0):.2e}, {r.get('intersection_size', 0)} genes: {gene_list_str})"
                    )

            if kegg_results:
                context_parts.append(f"KEGG Pathway Enrichment ({len(kegg_results)} pathways):")
                for r in kegg_results[:8]:
                    genes_in = r.get("genes", [])
                    gene_list_str = ", ".join(genes_in[:5]) if isinstance(genes_in, list) else str(genes_in)[:60]
                    context_parts.append(
                        f"  - {r.get('term', 'N/A')} (p={r.get('pvalue', 0):.2e}, {r.get('intersection_size', 0)} genes: {gene_list_str})"
                    )

        if not context_parts:
            return jsonify({"error": "No analysis data available. Please run analysis first."}), 400

        user_message = "\n".join(context_parts)
        system_prompt = (
            "You are an expert bioinformatician. Based on the following gene analysis data, "
            "provide a comprehensive summary in English. The summary should include:\n"
            "1. Overview of the gene lists and their sizes\n"
            "2. Key intersection patterns and shared genes\n"
            "3. Biological interpretation of GO enrichment results\n"
            "4. Biological interpretation of KEGG pathway results\n"
            "5. Key insights and potential biological significance\n"
            "6. Suggestions for further analysis\n\n"
            "Be specific, cite gene names, p-values, and pathways. Write in a professional academic tone."
        )

        payload = {
            "model": DEEPSEEK_MODEL,
            "messages": [
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_message}
            ],
            "temperature": 0.7,
            "max_tokens": 20000
        }

        resp = requests.post(
            DEEPSEEK_API_URL,
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {DEEPSEEK_API_KEY}"
            },
            json=payload,
            timeout=60
        )

        if resp.status_code != 200:
            logger.error(f"DeepSeek API error: {resp.status_code} {resp.text[:200]}")
            return jsonify({"error": f"AI API error: {resp.status_code}"}), 500

        result = resp.json()
        summary = result.get("choices", [{}])[0].get("message", {}).get("content", "No summary generated.")

        return jsonify({"summary": summary})

    except requests.Timeout:
        return jsonify({"error": "AI request timed out. Please try again."}), 504
    except Exception as e:
        logger.error(f"AI summary failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


if __name__ == "__main__":
    # Disable debug in production, use gunicorn or other WSGI server
    app.run(host="0.0.0.0", port=8001, debug=True)
