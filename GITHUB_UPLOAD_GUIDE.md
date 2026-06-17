# Bioinformatics Platform GitHub Upload Guide

> This package is a sanitized code-only export of the 43 server services running on ports 8000-8006. It is intended for GitHub submission/review, not as a full production data dump.

## 1. Package contents

| Directory | Port | Project | Notes |
|---|---:|---|---|
| `animal_hormone_8000/` | 8000 | AHormoneDB / animal hormone database | Code only; BLAST DB, RAG embeddings, local models, media/uploads excluded |
| `venn_tool_8001/` | 8001 | Venn / single-cell visualization tool | Code only; single-cell datasets excluded |
| `customer_service_8002/` | 8002 | Bioinformatics customer service chatbot | Code only |
| `gene_reg_8003/` | 8003 | GeneReg transcription factor-target prediction | Code only; ENCODE peak data excluded |
| `multi_omics_8004/` | 8004 | Multi-omics analysis platform | Code only; runtime data/results excluded |
| `literature_search_8005/` | 8005 | Literature search / RAG helper | Code only |
| `gene_mapper_8006/` | 8006 | Gene mapper | Code only |

## 2. Sanitization performed

Before packaging, the export copy was cleaned as follows:

- Removed secret/runtime files: `.env`, `*.key`, `*.pem`, `*.p12`, `*.pfx`, `*.sqlite3`, `*.db`, `*.log`, `*.pid`.
- Removed backups and cache: `*.bak*`, `backups/`, `__pycache__/`, `venv/`, `node_modules/`.
- Removed large data/runtime directories: `blast/`, `local_models/`, `media/`, `uploads/`, `single_cell_data/`, `single_cell_data_mouse/`, `encode_peaks/`, `data/`, `results/`.
- Replaced hardcoded sensitive values with placeholders such as `YOUR_VALUE_HERE` and `YOUR_API_KEY_HERE`.
- Added `.env.example` files for each project.
- Added root and per-project `.gitignore` files.

Final verification status:

```text
Forbidden file scan: OK
Forbidden directory scan: OK
High-confidence secret scan: TOTAL_REMAINING_HIGH_CONFIDENCE=0
```

## 3. Important warning

Do **not** upload the old unsafe package:

```text
/root/bioinformatics_platform_github_20260618_UNSAFE_DO_NOT_UPLOAD.tar.gz
```

Use only the safe package:

```text
/root/bioinformatics_platform_github_20260618_SAFE.tar.gz
```

## 4. Create a GitHub repository and upload

### Option A: Upload from your local computer

Download the safe package from the server:

```bash
scp root@43.99.62.219:/root/bioinformatics_platform_github_20260618_SAFE.tar.gz .
```

Extract it:

```bash
mkdir bioinformatics_platform
cd bioinformatics_platform
tar -xzf ../bioinformatics_platform_github_20260618_SAFE.tar.gz
```

Initialize Git and push:

```bash
git init
git add .
git commit -m "Initial sanitized bioinformatics platform export"
git branch -M main
git remote add origin https://github.com/<your-username>/<your-repo>.git
git push -u origin main
```

### Option B: Upload from a server with GitHub CLI

```bash
mkdir bioinformatics_platform
cd bioinformatics_platform
tar -xzf /root/bioinformatics_platform_github_20260618_SAFE.tar.gz
git init
git add .
git commit -m "Initial sanitized bioinformatics platform export"
gh repo create <your-repo> --private --source=. --remote=origin --push
```

Use `--public` only if you are sure the repository can be public.

## 5. Configure secrets after cloning

For each service, copy `.env.example` to `.env` and fill real values locally/server-side only:

```bash
cp .env.example .env
```

Never commit `.env`.

Typical required variables include:

- `SECRET_KEY`
- `DATABASE_NAME`
- `DATABASE_USER`
- `DATABASE_PASSWORD`
- `DATABASE_HOST`
- `DEEPSEEK_API_KEY`
- `SERPER_API_KEY`
- `OPENAI_API_KEY` if used

## 6. Restore large data files if deploying

The GitHub package does not include heavy data. Production deployment requires restoring these separately:

- AHormoneDB: BLAST databases, RAG embeddings/text/metadata, local embedding model if used.
- Venn tool: single-cell datasets such as `.h5ad`/`.zarr` files.
- GeneReg: ENCODE ChIP-seq peak files under `encode_peaks/`.
- Multi-omics: local data/results directories as needed.

Recommended practice: store these in object storage, institutional storage, or a separate private data archive, not GitHub.

## 7. Suggested repository hygiene

Before every public push, run:

```bash
git status
find . -type f \( -name '.env' -o -name '*.key' -o -name '*.pem' -o -name '*.sqlite3' -o -name '*.db' \)
git grep -n -i -E "(SECRET_KEY|api[_-]?key|token|password|passwd|private key)" || true
```

If possible, install and run a secret scanner such as `gitleaks` or `trufflehog` before making the repository public.

## 8. Current safe package metadata

```text
File: bioinformatics_platform_github_20260618_SAFE.tar.gz
Extracted size: about 31 MB
SHA256: 768680eec7ef8ce8dedb897f344cdc97e260c1f0262e49539dd0b14c79db1116
```


## 9. License and citation

This repository uses the MIT License. See `LICENSE` for details.

If you use or reference this codebase, please cite the repository using `CITATION.cff`.
