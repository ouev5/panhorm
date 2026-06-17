#!/usr/bin/env python3
"""
Gene Mapper - Cross-species gene ortholog mapping tool
Maps genes from any species to human orthologs using NCBI HomoloGene & MyGene.info APIs
"""

from flask import Flask, render_template, request, jsonify
import requests
import logging
import sys
import os
import time
import pymysql

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 10 * 1024 * 1024

# Species info: NCBI taxid -> name
SPECIES_MAP = {
    9606: "Human",
    10090: "Mouse",
    10116: "Rat",
    9913: "Bovine",
    9823: "Pig",
    9940: "Sheep",
    9031: "Chicken",
    9615: "Dog",
    9986: "Rabbit",
    9685: "Cat",
    8364: "African clawed frog",
    9796: "Horse",
    9544: "Rhesus macaque",
    10141: "Guinea pig",
    9541: "Crab-eating macaque",
    7227: "Fruit fly",
    8049: "Rainbow trout",
    7955: "Zebrafish",
    9598: "Chimpanzee",
    9601: "Sumatran orangutan",
    8432: "American bullfrog",
    9925: "Goat",
    10036: "Golden hamster",
    7994: "Common carp",
    9103: "Wild turkey",
    7998: "Goldfish",
    8258: "Gilthead sea bream",
    8030: "Atlantic salmon",
    9483: "White-tufted-ear marmoset",
    8187: "Bastard halibut",
}

# Reverse map: name -> taxid
SPECIES_NAME_TO_ID = {v.lower(): k for k, v in SPECIES_MAP.items()}
SPECIES_NAME_TO_ID["human"] = 9606
SPECIES_NAME_TO_ID["mouse"] = 10090
SPECIES_NAME_TO_ID["rat"] = 10116
SPECIES_NAME_TO_ID["bovine"] = 9913
SPECIES_NAME_TO_ID["pig"] = 9823
SPECIES_NAME_TO_ID["sheep"] = 9940
SPECIES_NAME_TO_ID["chicken"] = 9031
SPECIES_NAME_TO_ID["dog"] = 9615
SPECIES_NAME_TO_ID["rabbit"] = 9986
SPECIES_NAME_TO_ID["cat"] = 9685
SPECIES_NAME_TO_ID["horse"] = 9796
SPECIES_NAME_TO_ID["zebrafish"] = 7955
SPECIES_NAME_TO_ID["fruit fly"] = 7227
SPECIES_NAME_TO_ID["drosophila"] = 7227
SPECIES_NAME_TO_ID["c. elegans"] = 6239
SPECIES_NAME_TO_ID["nematode"] = 6239

# Popular species for dropdown
POPULAR_SPECIES = [9606, 10090, 10116, 7955, 7227, 9913, 9823, 9940, 9031, 9615, 9986, 9685]


def get_db_connection():
    try:
        conn = pymysql.connect(
            host='127.0.0.1', user='hhh', password='YOUR_VALUE_HERE',
            database='animal_hormone', charset='utf8mb4',
            cursorclass=pymysql.cursors.DictCursor
        )
        return conn
    except Exception as e:
        logger.error(f"DB connection failed: {e}")
        return None


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/species", methods=["GET"])
def get_species():
    """Return available species list"""
    result = []
    for taxid in POPULAR_SPECIES:
        if taxid in SPECIES_MAP:
            result.append({"taxid": taxid, "name": SPECIES_MAP[taxid]})
    # Add all others
    for taxid, name in sorted(SPECIES_MAP.items(), key=lambda x: x[1]):
        if taxid not in POPULAR_SPECIES:
            result.append({"taxid": taxid, "name": name})
    return jsonify({"species": result})


@app.route("/api/available_species", methods=["GET"])
def available_species():
    """Get species from database that have gene data"""
    try:
        conn = get_db_connection()
        if not conn:
            return jsonify({"species": []})
        cur = conn.cursor()
        cur.execute("""
            SELECT DISTINCT organism, COUNT(*) as cnt 
            FROM hormone_data 
            WHERE organism IS NOT NULL AND organism != '' AND related_genes IS NOT NULL AND related_genes != ''
            GROUP BY organism ORDER BY cnt DESC
        """)
        db_species = []
        for row in cur.fetchall():
            db_species.append({"name": row["organism"], "count": row["cnt"]})
        cur.close()
        conn.close()
        return jsonify({"species": db_species})
    except Exception as e:
        return jsonify({"species": [], "error": str(e)})


@app.route("/api/map_genes", methods=["POST"])
def map_genes():
    """Map genes from source species to human orthologs"""
    try:
        data = request.get_json(force=True)
        genes = data.get("genes", [])
        source_species = data.get("source_species", "")
        target_taxid = data.get("target_taxid", 9606)  # default: human

        if not genes:
            return jsonify({"error": "No genes provided", "results": []}), 400

        # Clean gene list
        genes = [g.strip().upper() for g in genes if g.strip()]
        if not genes:
            return jsonify({"error": "No valid genes", "results": []}), 400

        # Determine source taxid
        source_taxid = None
        if isinstance(source_species, int):
            source_taxid = source_species
        elif isinstance(source_species, str):
            source_taxid = SPECIES_NAME_TO_ID.get(source_species.lower())
            if not source_taxid:
                # Try as string number
                try:
                    source_taxid = int(source_species)
                except ValueError:
                    pass

        # Method 1: Use MyGene.info API for batch ortholog mapping
        results = batch_map_mygene(genes, source_taxid, target_taxid)

        return jsonify({
            "success": True,
            "count": len(results),
            "results": results,
            "source_species": source_species,
            "target_taxid": target_taxid,
            "target_name": SPECIES_MAP.get(target_taxid, "Unknown")
        })

    except Exception as e:
        logger.error(f"Gene mapping failed: {e}")
        import traceback
        traceback.print_exc()
        return jsonify({"error": str(e), "results": []}), 500


def batch_map_mygene(genes, source_taxid, target_taxid):
    """Use MyGene.info API to batch map genes to orthologs"""
    results = []
    
    # Batch query MyGene.info (max 1000 per request)
    batch_size = 50
    for i in range(0, len(genes), batch_size):
        batch = genes[i:i + batch_size]
        
        # Build query
        query_str = ",".join(batch)
        
        params = {
            "q": query_str,
            "scopes": "symbol,alias,name",
            "fields": "symbol,name,taxid,homologene",
            "size": len(batch) * 2,
        }
        
        if source_taxid:
            params["species"] = source_taxid
        
        try:
            resp = requests.post(
                "https://mygene.info/v3/query",
                data=params,
                timeout=30
            )
            
            if resp.status_code != 200:
                logger.error(f"MyGene API error: {resp.status_code}")
                for g in batch:
                    results.append({"query": g, "mapped_gene": None, "status": "api_error"})
                continue
            
            data = resp.json()
            
            # Build lookup: query -> best match (case-insensitive)
            gene_results = {}
            batch_upper = {g.upper(): g for g in batch}  # uppercase -> original

            for hit in data:
                if not isinstance(hit, dict):
                    continue
                
                symbol = hit.get("symbol", "").upper()
                taxid = hit.get("taxid")
                
                # Find original query gene name
                orig_query = batch_upper.get(symbol)
                if not orig_query:
                    continue

                # Check if this gene matches target species directly
                if taxid == target_taxid:
                    if orig_query not in gene_results:
                        gene_results[orig_query] = {
                            "query": orig_query,
                            "mapped_gene": hit.get("symbol", orig_query),
                            "mapped_name": hit.get("name", ""),
                            "status": "mapped",
                            "method": "direct_match",
                            "source_taxid": taxid
                        }
                    continue

                # Check HomoloGene for ortholog mapping
                homologene = hit.get("homologene")
                if not homologene or not isinstance(homologene, dict):
                    continue

                hg_genes = homologene.get("genes", [])
                target_entrez = None
                for hg in hg_genes:
                    if isinstance(hg, list) and len(hg) >= 2:
                        if hg[0] == target_taxid:
                            target_entrez = hg[1]
                            break
                    elif isinstance(hg, dict):
                        if hg.get("taxid") == target_taxid:
                            target_entrez = hg.get("geneid", hg.get("entrezgene"))
                            break

                if target_entrez and orig_query not in gene_results:
                    gene_results[orig_query] = {
                        "query": orig_query,
                        "mapped_gene": None,
                        "mapped_entrez": target_entrez,
                        "mapped_name": "",
                        "status": "mapped",
                        "method": "homologene",
                        "source_taxid": taxid
                    }

            # Fill unmapped genes
            for g in batch:
                if g not in gene_results:
                    gene_results[g] = {
                        "query": g,
                        "mapped_gene": None,
                        "status": "not_found"
                    }
            
            results.extend(gene_results.values())
            
        except requests.Timeout:
            for g in batch:
                results.append({"query": g, "mapped_gene": None, "status": "timeout"})
        except Exception as e:
            logger.error(f"MyGene batch error: {e}")
            for g in batch:
                results.append({"query": g, "mapped_gene": None, "status": "error"})
    
    # Enrich: convert mapped_entrez IDs to gene symbols via MyGene.info
    need_symbol = {str(r["mapped_entrez"]): i for i, r in enumerate(results)
                   if r.get("status") == "mapped" and r.get("mapped_entrez") and not r.get("mapped_gene")}
    
    if need_symbol:
        try:
            entrez_str = ",".join(list(need_symbol.keys())[:200])
            sym_resp = requests.post(
                "https://mygene.info/v3/query",
                data={
                    "q": entrez_str,
                    "scopes": "entrezgene",
                    "fields": "symbol,name,taxid",
                    "species": target_taxid,
                    "size": len(need_symbol) * 2
                },
                timeout=20
            )
            if sym_resp.status_code == 200:
                for hit in sym_resp.json():
                    if isinstance(hit, dict) and hit.get("symbol"):
                        eid = str(hit.get("entrezgene", hit.get("_id", "")))
                        # Try to match
                        for eid_key, idx in list(need_symbol.items()):
                            if eid_key == eid or eid_key == str(hit.get("_id", "")):
                                results[idx]["mapped_gene"] = hit["symbol"]
                                results[idx]["mapped_name"] = hit.get("name", "")
                                break
        except Exception as e:
            logger.error(f"Symbol lookup failed: {e}")

    # Also enrich direct matches that might lack names
    need_name = [i for i, r in enumerate(results)
                 if r.get("status") == "mapped" and r.get("mapped_gene") and not r.get("mapped_name")]
    
    if need_name:
        try:
            syms = [results[i]["mapped_gene"] for i in need_name[:200]]
            sym_str = ",".join(syms)
            nm_resp = requests.post(
                "https://mygene.info/v3/query",
                data={
                    "q": sym_str,
                    "scopes": "symbol",
                    "fields": "symbol,name",
                    "species": target_taxid,
                    "size": len(syms) * 2
                },
                timeout=15
            )
            if nm_resp.status_code == 200:
                name_map = {}
                for hit in nm_resp.json():
                    if isinstance(hit, dict) and hit.get("symbol"):
                        name_map[hit["symbol"].upper()] = hit.get("name", "")
                for idx in need_name:
                    r = results[idx]
                    if not r.get("mapped_name"):
                        r["mapped_name"] = name_map.get(r["mapped_gene"].upper(), "")
        except Exception:
            pass

    # Fallback: if mapped_gene still None, use entrez id as string
    for r in results:
        if r.get("status") == "mapped" and not r.get("mapped_gene") and r.get("mapped_entrez"):
            r["mapped_gene"] = str(r["mapped_entrez"])

    return results


@app.route("/api/map_from_db", methods=["POST"])
def map_from_db():
    """Map genes from a specific species in our database to target species"""
    try:
        data = request.get_json(force=True)
        source_organism = data.get("source_organism", "")
        target_taxid = data.get("target_taxid", 9606)
        gene_list = data.get("genes", [])

        if not source_organism:
            return jsonify({"error": "Source organism required", "results": []}), 400

        conn = get_db_connection()
        if not conn:
            return jsonify({"error": "Database connection failed", "results": []}), 500

        cur = conn.cursor()
        
        # Get all genes for this organism
        if gene_list:
            placeholders = ",".join(["%s"] * len(gene_list))
            cur.execute(f"""
                SELECT DISTINCT related_genes FROM hormone_data 
                WHERE organism = %s AND related_genes != '' AND related_genes IS NOT NULL
                AND related_genes IN ({placeholders})
            """, [source_organism] + gene_list)
        else:
            cur.execute("""
                SELECT DISTINCT related_genes FROM hormone_data 
                WHERE organism = %s AND related_genes != '' AND related_genes IS NOT NULL
                LIMIT 5000
            """, (source_organism,))

        all_genes = set()
        for row in cur.fetchall():
            for g in str(row["related_genes"]).split(","):
                g = g.strip()
                if g:
                    all_genes.add(g.upper())

        cur.close()
        conn.close()

        if not all_genes:
            return jsonify({"error": "No genes found", "results": []}), 400

        # Map to target
        gene_list_sorted = sorted(all_genes)
        source_taxid = SPECIES_NAME_TO_ID.get(source_organism.lower())
        results = batch_map_mygene(gene_list_sorted, source_taxid, target_taxid)

        return jsonify({
            "success": True,
            "source_organism": source_organism,
            "total_genes": len(gene_list_sorted),
            "mapped_count": sum(1 for r in results if r.get("status") == "mapped"),
            "results": results,
            "target_name": SPECIES_MAP.get(target_taxid, "Unknown")
        })

    except Exception as e:
        logger.error(f"DB gene mapping failed: {e}")
        import traceback
        traceback.print_exc()
        return jsonify({"error": str(e), "results": []}), 500


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8006, debug=True)
