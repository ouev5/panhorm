#!/usr/bin/env python3
"""
Multi-Omics Analysis Platform - Production Grade
Self-contained real analysis using Python scientific libraries
"""
from flask import Flask, render_template, request, jsonify
import os, json, uuid, traceback
import pandas as pd
import numpy as np
from datetime import datetime
from pathlib import Path

app = Flask(__name__)
app.config['SECRET_KEY'] = 'multi-omics-platform-2024'
app.config['UPLOAD_FOLDER'] = '/www/wwwroot/multi-omics/data'
app.config['OUTPUT_FOLDER'] = '/www/wwwroot/multi-omics/output'
app.config['MAX_CONTENT_LENGTH'] = 500 * 1024 * 1024

DEEPSEEK_API_KEY = "YOUR_VALUE_HERE"
DEEPSEEK_API_URL = "https://api.deepseek.com/v1/chat/completions"

# --- Pipeline definitions ---
PIPELINES = {
    'rnaseq_diffexpr': {
        'name': 'RNA-seq Differential Expression',
        'icon': 'fa-dna',
        'description': 'Identify differentially expressed genes from count/expression matrix',
        'params': ['pvalue_cutoff', 'log2fc_cutoff', 'method'],
        'accepts': 'CSV/TSV with genes × samples (first column = gene names, header row = sample names, last row of header can have group labels)',
    },
    'scrna_cluster': {
        'name': 'Single-cell RNA-seq Clustering',
        'icon': 'fa-circle-nodes',
        'description': 'Cell type identification and UMAP visualization',
        'params': ['n_pcs', 'resolution'],
        'accepts': 'CSV/TSV with genes × cells count matrix',
    },
    'proteomics_quant': {
        'name': 'Proteomics Quantification',
        'icon': 'fa-atom',
        'description': 'Protein expression analysis and normalization',
        'params': ['normalization', 'missing_value'],
        'accepts': 'CSV/TSV with proteins × samples intensity matrix',
    },
    'volcano': {
        'name': 'Volcano Plot',
        'icon': 'fa-fire',
        'description': 'Generate volcano plot from DE results or expression data',
        'params': ['pvalue_cutoff', 'log2fc_cutoff'],
        'accepts': 'CSV with columns: gene, log2FC, pvalue  OR  expression matrix with group info',
    },
    'correlation': {
        'name': 'Correlation Analysis',
        'icon': 'fa-link',
        'description': 'Gene-gene or feature correlation heatmap',
        'params': ['method'],
        'accepts': 'CSV/TSV expression matrix',
    },
    'deseq2': {
        'name': 'DESeq2 Differential Expression (R)',
        'icon': 'fa-bar-chart',
        'description': 'Publication-grade DE analysis using R DESeq2 with volcano/MA/heatmap plots',
        'params': ['pvalue_cutoff', 'log2fc_cutoff', 'control_prefix', 'treat_prefix'],
        'accepts': 'CSV count matrix (genes x samples, first col = gene names)',
    },
}

# --- Task storage ---
tasks = {}

# ========================
# DeepSeek Chat
# ========================
def call_deepseek(prompt, system_prompt=""):
    import requests
    try:
        headers = {"Authorization": f"Bearer {DEEPSEEK_API_KEY}", "Content-Type": "application/json"}
        messages = []
        if system_prompt:
            messages.append({"role": "system", "content": system_prompt})
        messages.append({"role": "user", "content": prompt})
        resp = requests.post(DEEPSEEK_API_URL, headers=headers,
                             json={"model": "deepseek-chat", "messages": messages,
                                   "temperature": 0.7, "max_tokens": 2000}, timeout=60)
        resp.raise_for_status()
        return resp.json()["choices"][0]["message"]["content"]
    except Exception as e:
        return f"Error: {str(e)}"

# ========================
# File parsing
# ========================
def parse_uploaded_file(file_path):
    ext = Path(file_path).suffix.lower()
    if ext == '.csv':
        return pd.read_csv(file_path)
    elif ext in ['.tsv', '.txt']:
        return pd.read_csv(file_path, sep='\t')
    elif ext in ['.xlsx', '.xls']:
        return pd.read_excel(file_path)
    else:
        raise ValueError(f"Unsupported format: {ext}")

# ========================
# Real analysis functions
# ========================
def do_volcano(df, pval_cutoff=0.05, log2fc_cutoff=1.0):
    """Real volcano plot computation.
    Accepts:
      - Direct: columns gene, log2FC, pvalue
      - Expression matrix: genes × samples, with optional group row
    """
    # Direct format
    if 'log2FC' in df.columns and 'pvalue' in df.columns:
        gene_col = df['gene'] if 'gene' in df.columns else pd.Series([f'Gene_{i}' for i in range(len(df))])
        result = pd.DataFrame({
            'gene': gene_col,
            'log2FC': df['log2FC'].astype(float),
            'pvalue': df['pvalue'].astype(float),
        })
    else:
        # Expression matrix: try to compute DE
        # Detect gene column
        gene_col = None
        for col in ['gene', 'Gene', 'gene_id', 'gene_name', 'Symbol', 'symbol']:
            if col in df.columns:
                gene_col = df[col]
                df = df.drop(columns=[col])
                break

        numeric_df = df.select_dtypes(include=[np.number])
        if numeric_df.empty:
            raise ValueError("No numeric columns found for analysis")

        ncols = len(numeric_df.columns)
        if ncols < 2:
            raise ValueError("Need at least 2 samples for differential expression")

        # Split: first half = control, second half = treatment
        mid = ncols // 2
        if mid == 0:
            mid = 1
        ctrl_cols = numeric_df.columns[:mid]
        treat_cols = numeric_df.columns[mid:]

        from scipy import stats
        log2fc_list = []
        pval_list = []
        for i in range(len(numeric_df)):
            ctrl_vals = numeric_df.iloc[i][ctrl_cols].values.astype(float)
            treat_vals = numeric_df.iloc[i][treat_cols].values.astype(float)

            # Add small pseudocount
            ctrl_vals = ctrl_vals + 1
            treat_vals = treat_vals + 1

            fc = np.mean(treat_vals) / np.mean(ctrl_vals)
            log2fc = np.log2(fc) if fc > 0 else 0

            if len(ctrl_vals) >= 2 and len(treat_vals) >= 2:
                _, pval = stats.ttest_ind(ctrl_vals, treat_vals)
                pval = pval if not np.isnan(pval) else 1.0
            else:
                pval = 1.0

            log2fc_list.append(log2fc)
            pval_list.append(pval)

        genes = gene_col.tolist() if gene_col is not None else [f'Gene_{i}' for i in range(len(numeric_df))]
        result = pd.DataFrame({'gene': genes, 'log2FC': log2fc_list, 'pvalue': pval_list})

    result['-log10(p)'] = -np.log10(result['pvalue'].clip(lower=1e-300))
    result['significant'] = ((abs(result['log2FC']) > log2fc_cutoff) & (result['pvalue'] < pval_cutoff)).astype(int)

    # Summary stats
    n_total = len(result)
    n_up = int(((result['log2FC'] > log2fc_cutoff) & (result['pvalue'] < pval_cutoff)).sum())
    n_down = int(((result['log2FC'] < -log2fc_cutoff) & (result['pvalue'] < pval_cutoff)).sum())
    n_sig = n_up + n_down

    return result.to_dict(orient='records'), {
        'total_genes': n_total, 'significant': n_sig,
        'up_regulated': n_up, 'down_regulated': n_down,
        'pval_cutoff': pval_cutoff, 'log2fc_cutoff': log2fc_cutoff
    }


def do_umap(df, n_pcs=30, resolution=0.5):
    """Real UMAP computation using Scanpy."""
    import scanpy as sc

    # Detect gene column
    gene_col = None
    for col in ['gene', 'Gene', 'gene_id', 'gene_name']:
        if col in df.columns:
            gene_col = df[col]
            df = df.drop(columns=[col])
            break

    numeric_df = df.select_dtypes(include=[np.number])
    if numeric_df.empty:
        raise ValueError("No numeric columns found")

    # Build AnnData: rows = genes, cols = cells → transpose
    adata = sc.AnnData(numeric_df.T.astype(float))
    if gene_col is not None:
        adata.var_names = gene_col.values

    # Adaptive filtering (lower thresholds for small datasets)
    n_cells = adata.n_obs
    min_genes_per_cell = max(1, min(200, n_cells))
    sc.pp.filter_cells(adata, min_genes=min_genes_per_cell)
    sc.pp.filter_genes(adata, min_cells=max(1, min(3, adata.n_obs // 10)))

    if adata.n_obs < 3:
        raise ValueError(f"Too few cells after filtering ({adata.n_obs}). Need at least 3 cells.")
    if adata.n_vars < 3:
        raise ValueError(f"Too few genes after filtering ({adata.n_vars}). Need at least 3 genes.")

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG - use 'seurat' flavor (not v3) for small datasets
    n_top = min(2000, adata.n_vars - 1)
    if n_top > 10:
        flavor = 'seurat' if n_cells < 100 else 'seurat_v3'
        try:
            sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor=flavor, subset=True)
        except Exception:
            pass  # Skip HVG if it fails

    # PCA
    n_comps = min(n_pcs, min(adata.n_vars, adata.n_obs) - 1)
    if n_comps < 2:
        n_comps = 2
    sc.tl.pca(adata, n_comps=n_comps)

    # Neighbors + UMAP
    n_nb = min(15, max(2, adata.n_obs - 1))
    sc.pp.neighbors(adata, n_neighbors=n_nb, n_pcs=min(n_comps, adata.n_obs - 1))
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution)

    result = pd.DataFrame({
        'UMAP1': adata.obsm['X_umap'][:, 0],
        'UMAP2': adata.obsm['X_umap'][:, 1],
        'cluster': adata.obs['leiden'].astype(str).values,
        'cell': adata.obs_names,
    })

    n_clusters = result['cluster'].nunique()
    return result.to_dict(orient='records'), {
        'n_cells': len(result), 'n_clusters': n_clusters,
        'n_genes_used': adata.n_vars
    }


def do_proteomics(df, normalization='median', missing_threshold=0.5):
    """Proteomics quantification analysis."""
    # Detect protein column
    protein_col = None
    for col in ['Protein', 'protein', 'ProteinID', 'Accession', 'Gene', 'gene']:
        if col in df.columns:
            protein_col = df[col]
            df = df.drop(columns=[col])
            break

    numeric_df = df.select_dtypes(include=[np.number])
    if numeric_df.empty:
        raise ValueError("No numeric columns found")

    # Filter by missing value threshold
    valid_mask = numeric_df.isnull().mean(axis=1) < missing_threshold
    numeric_df = numeric_df[valid_mask]
    proteins = protein_col[valid_mask] if protein_col is not None else pd.Series([f'Protein_{i}' for i in range(len(numeric_df))])

    # Impute missing values with minimum/2
    numeric_df = numeric_df.fillna(numeric_df.min().min() / 2)

    # Normalization
    if normalization == 'median':
        medians = numeric_df.median()
        numeric_df = numeric_df.divide(medians) * numeric_df.median().median()
    elif normalization == 'quantile':
        from scipy.stats import rankdata
        rank_mean = numeric_df.stack().groupby(rankdata(numeric_df.values.flatten())).mean()
        for col in numeric_df.columns:
            ranks = rankdata(numeric_df[col].values)
            numeric_df[col] = [rank_mean[r] for r in ranks]

    # Log2 transform
    log2_df = np.log2(numeric_df.clip(lower=0.01))

    # Stats
    stats = {
        'n_proteins': len(log2_df),
        'n_samples': len(log2_df.columns),
        'normalization': normalization,
        'missing_filtered': int((~valid_mask).sum()),
    }

    # Return top variable proteins for heatmap
    variances = log2_df.var(axis=1)
    top_idx = variances.nlargest(min(100, len(variances))).index
    heatmap_data = log2_df.loc[top_idx]

    heatmap_result = []
    for idx, row in heatmap_data.iterrows():
        for col in row.index:
            heatmap_result.append({
                'protein': proteins.loc[idx] if idx in proteins.index else f'Protein_{idx}',
                'sample': col,
                'value': round(float(row[col]), 3)
            })

    return heatmap_result, stats


def do_correlation(df, method='pearson', top_n=50):
    """Gene-gene correlation analysis."""
    # Detect gene column
    for col in ['gene', 'Gene', 'gene_id', 'gene_name']:
        if col in df.columns:
            df = df.drop(columns=[col])
            break

    numeric_df = df.select_dtypes(include=[np.number])
    if numeric_df.empty:
        raise ValueError("No numeric columns found")

    # Take top N most variable genes (rows)
    variances = numeric_df.var(axis=1)
    top_idx = variances.nlargest(min(top_n, len(variances))).index
    sub = numeric_df.loc[top_idx]

    # Correlation matrix
    if method == 'spearman':
        corr = sub.T.corr(method='spearman')
    else:
        corr = sub.T.corr(method='pearson')

    # Build result as list of {gene1, gene2, corr}
    genes = corr.index.tolist()
    result = []
    for i, g1 in enumerate(genes):
        for j, g2 in enumerate(genes):
            if i < j:  # upper triangle only
                result.append({'gene1': str(g1), 'gene2': str(g2), 'corr': round(float(corr.iloc[i, j]), 4)})

    return result, {
        'n_genes': len(genes), 'method': method,
        'mean_abs_corr': round(float(corr.abs().mean().mean()), 4)
    }


# ========================
# Routes
# ========================
@app.route('/')
def index():
    return render_template('index.html', pipelines=PIPELINES)

@app.route('/api/chat', methods=['POST'])
def chat():
    try:
        message = request.json.get('message', '')
        system_prompt = """You are an AI assistant for a multi-omics bioinformatics analysis platform.
Respond in English with proper Markdown formatting (## headers, bullet lists, **bold**, code blocks).
Help researchers with RNA-seq, scRNA-seq, proteomics, ChIP-seq, and ATAC-seq analysis.
Provide detailed, professional guidance with clear structure."""
        response = call_deepseek(message, system_prompt)
        return jsonify({'response': response})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/upload', methods=['POST'])
def upload_file():
    try:
        if 'file' not in request.files:
            return jsonify({'error': 'No file provided'}), 400
        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'Empty filename'}), 400
        filename = f"{uuid.uuid4().hex}_{file.filename}"
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
        file.save(filepath)

        # Quick preview
        try:
            df = parse_uploaded_file(filepath)
            preview = {
                'rows': len(df), 'cols': len(df.columns),
                'columns': list(df.columns[:20]),
                'head': df.head(3).to_dict(orient='records')
            }
        except:
            preview = {'rows': '?', 'cols': '?', 'columns': [], 'head': []}

        return jsonify({'success': True, 'filename': filename, 'filepath': filepath, 'preview': preview})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/pipelines')
def get_pipelines():
    return jsonify(PIPELINES)

@app.route('/api/run_analysis', methods=['POST'])
def run_analysis():
    """Run real analysis on uploaded data."""
    try:
        data = request.json
        task_id = str(uuid.uuid4())[:8]
        pipeline = data.get('pipeline', 'rnaseq_diffexpr')
        filepath = data.get('filepath', '')
        params = {}

        # Parse advanced params
        if data.get('advanced'):
            try:
                params = json.loads(data['advanced'])
            except:
                pass

        tasks[task_id] = {'status': 'running', 'pipeline': pipeline, 'started': datetime.now().isoformat()}
        logs = [f"[{datetime.now().strftime('%H:%M:%S')}] Task {task_id} started",
                f"[{datetime.now().strftime('%H:%M:%S')}] Pipeline: {pipeline}"]

        if not filepath or not os.path.exists(filepath):
            logs.append(f"[{datetime.now().strftime('%H:%M:%S')}] ERROR: No data file uploaded. Please upload a file first.")
            tasks[task_id] = {'status': 'failed', 'pipeline': pipeline, 'logs': logs}
            return jsonify({'task_id': task_id, 'status': 'failed', 'logs': logs})

        # Parse file
        logs.append(f"[{datetime.now().strftime('%H:%M:%S')}] Parsing uploaded file...")
        df = parse_uploaded_file(filepath)
        logs.append(f"[{datetime.now().strftime('%H:%M:%S')}] Loaded: {len(df)} rows × {len(df.columns)} columns")

        # Run analysis
        logs.append(f"[{datetime.now().strftime('%H:%M:%S')}] Running analysis...")
        output_dir = os.path.join(app.config['OUTPUT_FOLDER'], task_id)
        os.makedirs(output_dir, exist_ok=True)

        result_data = None
        stats = {}

        if pipeline == 'rnaseq_diffexpr' or pipeline == 'volcano':
            pval_cutoff = params.get('pvalue_cutoff', 0.05)
            fc_cutoff = params.get('log2fc_cutoff', 1.0)
            result_data, stats = do_volcano(df, pval_cutoff, fc_cutoff)
            # Save results
            pd.DataFrame(result_data).to_csv(os.path.join(output_dir, 'volcano_results.csv'), index=False)
            logs.append(f"[{datetime.now().strftime('%H:%M:%S')}] DE analysis: {stats['significant']} significant genes ({stats['up_regulated']} up, {stats['down_regulated']} down)")

        elif pipeline == 'scrna_cluster':
            n_pcs = params.get('n_pcs', 30)
            resolution = params.get('resolution', 0.5)
            result_data, stats = do_umap(df, n_pcs, resolution)
            pd.DataFrame(result_data).to_csv(os.path.join(output_dir, 'umap_results.csv'), index=False)
            logs.append(f"[{datetime.now().strftime('%H:%M:%S')}] scRNA-seq: {stats['n_cells']} cells, {stats['n_clusters']} clusters identified")

        elif pipeline == 'proteomics_quant':
            norm = params.get('normalization', 'median')
            result_data, stats = do_proteomics(df, norm)
            pd.DataFrame(result_data).to_csv(os.path.join(output_dir, 'proteomics_results.csv'), index=False)
            logs.append(f"[{datetime.now().strftime('%H:%M:%S')}] Proteomics: {stats['n_proteins']} proteins quantified across {stats['n_samples']} samples")

        elif pipeline == 'correlation':
            method = params.get('method', 'pearson')
            top_n = params.get('top_n', 50)
            result_data, stats = do_correlation(df, method, top_n)
            pd.DataFrame(result_data).to_csv(os.path.join(output_dir, 'correlation_results.csv'), index=False)
            logs.append(f"[{datetime.now().strftime('%H:%M:%S')}] Correlation: {stats['n_genes']} genes analyzed, mean |r| = {stats['mean_abs_corr']}")

        else:
            # Default: try volcano
            result_data, stats = do_volcano(df)
            pd.DataFrame(result_data).to_csv(os.path.join(output_dir, 'results.csv'), index=False)
            logs.append(f"[{datetime.now().strftime('%H:%M:%S')}] Default analysis completed")

        logs.append(f"[{datetime.now().strftime('%H:%M:%S')}] Results saved to {output_dir}")
        logs.append(f"[{datetime.now().strftime('%H:%M:%S')}] ✅ All steps completed successfully!")

        tasks[task_id] = {'status': 'completed', 'pipeline': pipeline, 'logs': logs, 'stats': stats}
        return jsonify({'task_id': task_id, 'status': 'completed', 'logs': logs,
                        'stats': stats, 'data': result_data, 'output_dir': output_dir})

    except Exception as e:
        traceback.print_exc()
        logs.append(f"[{datetime.now().strftime('%H:%M:%S')}] ❌ Error: {str(e)}")
        return jsonify({'task_id': task_id, 'status': 'failed', 'logs': logs, 'error': str(e)})


@app.route('/api/generate_volcano', methods=['POST'])
def generate_volcano():
    """Generate volcano plot data."""
    try:
        data = request.json or {}
        filepath = data.get('filepath')
        pval_cutoff = data.get('pvalue_cutoff', 0.05)
        fc_cutoff = data.get('log2fc_cutoff', 1.0)

        if not filepath or not os.path.exists(filepath):
            return jsonify({'error': 'No file uploaded. Please upload a CSV/TSV file with gene expression data.', 'needs_file': True}), 400

        df = parse_uploaded_file(filepath)
        result, stats = do_volcano(df, pval_cutoff, fc_cutoff)
        return jsonify({'data': result, 'stats': stats})

    except Exception as e:
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


@app.route('/api/generate_umap', methods=['POST'])
def generate_umap():
    """Generate UMAP data."""
    try:
        data = request.json or {}
        filepath = data.get('filepath')

        if not filepath or not os.path.exists(filepath):
            return jsonify({'error': 'No file uploaded. Please upload a CSV/TSV gene × cell count matrix.', 'needs_file': True}), 400

        df = parse_uploaded_file(filepath)
        result, stats = do_umap(df)
        return jsonify({'data': result, 'stats': stats})

    except Exception as e:
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


@app.route('/api/analyze_proteomics', methods=['POST'])
def analyze_proteomics():
    """Proteomics analysis."""
    try:
        data = request.json or {}
        filepath = data.get('filepath')
        normalization = data.get('normalization', 'median')

        if not filepath or not os.path.exists(filepath):
            return jsonify({'error': 'No file uploaded. Please upload a protein intensity matrix.', 'needs_file': True}), 400

        df = parse_uploaded_file(filepath)
        result, stats = do_proteomics(df, normalization)
        return jsonify({'data': result, 'stats': stats})

    except Exception as e:
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


@app.route('/api/analyze_correlation', methods=['POST'])
def analyze_correlation():
    """Correlation analysis."""
    try:
        data = request.json or {}
        filepath = data.get('filepath')
        method = data.get('method', 'pearson')
        top_n = data.get('top_n', 50)

        if not filepath or not os.path.exists(filepath):
            return jsonify({'error': 'No file uploaded.', 'needs_file': True}), 400

        df = parse_uploaded_file(filepath)
        result, stats = do_correlation(df, method, top_n)
        return jsonify({'data': result, 'stats': stats})

    except Exception as e:
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


# DESeq2 R Analysis (Publication-grade)
# ========================
@app.route('/api/run_deseq2', methods=['POST'])
def run_deseq2():
    """Run DESeq2 differential expression analysis using R."""
    try:
        import subprocess
        data = request.json or {}
        filepath = data.get('filepath')
        if not filepath or not os.path.exists(filepath):
            return jsonify({'error': 'No data file uploaded.'}), 400

        task_id = str(uuid.uuid4())[:8]
        output_dir = os.path.join(app.config['OUTPUT_FOLDER'], f'deseq2_{task_id}')
        os.makedirs(output_dir, exist_ok=True)

        # Parse options
        options = {}
        if data.get('pvalue_cutoff'): options['pvalue_cutoff'] = float(data['pvalue_cutoff'])
        if data.get('log2fc_cutoff'): options['log2fc_cutoff'] = float(data['log2fc_cutoff'])
        if data.get('control_prefix'): options['control_prefix'] = data['control_prefix']
        if data.get('treat_prefix'): options['treat_prefix'] = data['treat_prefix']
        options_json = json.dumps(options) if options else '{}'

        # Run R script
        r_script = os.path.join(os.path.dirname(__file__), 'scripts', 'run_deseq2.R')
        cmd = ['Rscript', r_script, filepath, output_dir, options_json]

        logs = [f"[{datetime.now().strftime('%H:%M:%S')}] Starting DESeq2 analysis (R)...",
                f"[{datetime.now().strftime('%H:%M:%S')}] Input: {filepath}"]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600, cwd=output_dir)

        # Parse R output
        r_logs = result.stdout.strip().split('\n') if result.stdout else []
        r_errors = result.stderr.strip().split('\n') if result.stderr else []
        for l in r_logs:
            if l.strip():
                logs.append(f"[{datetime.now().strftime('%H:%M:%S')}] {l.strip()}")

        if result.returncode != 0:
            for e in r_errors[-5:]:
                if 'Error' in e or 'error' in e:
                    logs.append(f"[{datetime.now().strftime('%H:%M:%S')}] R Error: {e.strip()}")
            logs.append(f"[{datetime.now().strftime('%H:%M:%S')}] ❌ DESeq2 analysis failed")
            return jsonify({'task_id': task_id, 'status': 'failed', 'logs': logs, 'error': 'R script failed'})

        # Read results
        results_csv = os.path.join(output_dir, 'deseq2_results.csv')
        summary_json = os.path.join(output_dir, 'summary.json')

        volcano_data = []
        stats = {}
        if os.path.exists(results_csv):
            res_df = pd.read_csv(results_csv)
            # For ECharts volcano: need log2FC, -log10(padj), gene, significant
            res_df['-log10(padj)'] = -np.log10(res_df['padj'].clip(lower=1e-300))
            volcano_data = res_df[['gene', 'log2FoldChange', 'pvalue', 'padj', '-log10(padj)', 'significant']].to_dict(orient='records')

        if os.path.exists(summary_json):
            with open(summary_json) as f:
                stats = json.load(f)

        # Image URLs
        images = {}
        for img_name in ['volcano_plot.png', 'ma_plot.png', 'heatmap.png']:
            img_path = os.path.join(output_dir, img_name)
            if os.path.exists(img_path):
                images[img_name] = f'/api/plot/{task_id}/{img_name}'

        logs.append(f"[{datetime.now().strftime('%H:%M:%S')}] ✅ DESeq2 analysis complete!")
        logs.append(f"[{datetime.now().strftime('%H:%M:%S')}] Results: {stats.get('total_genes', '?')} genes, {stats.get('significant', '?')} significant")

        return jsonify({
            'task_id': task_id,
            'status': 'completed',
            'logs': logs,
            'stats': stats,
            'data': volcano_data,
            'images': images,
            'output_dir': output_dir,
            'method': 'DESeq2 (R)'
        })

    except subprocess.TimeoutExpired:
        return jsonify({'error': 'DESeq2 analysis timed out (10 min limit)'}), 500
    except Exception as e:
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


@app.route('/api/plot/<task_id>/<filename>')
def serve_plot(task_id, filename):
    """Serve R-generated plot images."""
    from flask import send_from_directory
    output_dir = os.path.join(app.config['OUTPUT_FOLDER'], f'deseq2_{task_id}')
    if os.path.exists(os.path.join(output_dir, filename)):
        return send_from_directory(output_dir, filename, mimetype='image/png')
    return jsonify({'error': 'File not found'}), 404

if __name__ == '__main__':
    os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
    os.makedirs(app.config['OUTPUT_FOLDER'], exist_ok=True)
    app.run(host='0.0.0.0', port=8004, debug=False)
