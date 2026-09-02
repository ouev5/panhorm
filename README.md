# PanHorm Bioinformatics Platform

**PanHorm** is an integrated bioinformatics platform collection for hormone-related data mining, gene regulation analysis, multi-omics exploration, biomedical literature-assisted discovery, and interactive visualization.

This repository contains a **sanitized code-only export** of multiple web services originally deployed on a research server. It is designed for project demonstration, code review, academic sharing, and future reproducible development.

> **Security and data note**  
> Large datasets, API keys, local models, BLAST databases, ChIP-seq peak files, uploaded user files, runtime logs, and database files are intentionally excluded from this public repository for security, privacy, and storage reasons.

---

## Overview

PanHorm brings together several research-oriented web services used in hormone biology and bioinformatics workflows. The platform focuses on:

- Animal hormone data organization and hormone-gene-disease association exploration
- Transcription factor-target gene regulation prediction
- Literature-assisted biomedical evidence extraction
- Multi-omics workflow management and AI-assisted interpretation
- Venn, UpSet, enrichment, and single-cell visualization utilities
- Gene mapping and auxiliary analysis tools

The repository is arranged as a multi-service codebase. Each subdirectory corresponds to one original service/module.

## GeneReg v2.6 manuscript release

The `gene_reg_8003/` directory now includes the manuscript-matched GeneReg v2.6 snapshot and its frozen internal-benchmark materials. The release is identified as `AHormoneDB-GeneReg-v2.6-20260902`; see the module release README and manifest for checksums and audit files. Large runtime datasets and secrets remain excluded.

---

## Modules

| Directory | Original Port | Module | Description |
|---|---:|---|---|
| `animal_hormone_8000/` | 8000 | AHormoneDB / Animal Hormone Platform | Animal hormone database, hormone-gene-disease analysis, RAG-assisted Q&A, and hormone-related web interfaces |
| `venn_tool_8001/` | 8001 | Gene Insights / Venn Tool | Venn/UpSet visualization, enrichment analysis, and selected single-cell visualization utilities |
| `customer_service_8002/` | 8002 | Bioinformatics Assistant | AI-assisted bioinformatics chatbot/customer-service style module |
| `gene_reg_8003/` | 8003 | GeneReg v2.6 | Manuscript-matched transcription factor-target prediction with deterministic evidence scoring and frozen benchmark materials |
| `multi_omics_8004/` | 8004 | Multi-omics Platform | Multi-omics workflow interface and AI-assisted interpretation utilities |
| `literature_search_8005/` | 8005 | Literature Search | Biomedical literature search, evidence extraction, and RAG-style support utilities |
| `gene_mapper_8006/` | 8006 | Gene Mapper | Gene mapping and related analysis service |

---

## Repository Scope

This repository is intended to provide the **source code structure and implementation logic** of the platform. It is not a full production backup.

Included:

- Application source code
- Templates and static assets needed for code review
- Example environment files: `.env.example`
- Git ignore rules
- Upload and deployment notes
- Sanitized configuration examples

Excluded:

- `.env` files and real API keys
- Database files such as `*.db` and `*.sqlite3`
- User uploads and media files
- Runtime logs, PID files, and backup files
- Python virtual environments and dependency caches
- `node_modules/`
- BLAST databases
- Local embedding models
- RAG embeddings, metadata, and generated indexes
- ENCODE ChIP-seq peak files
- Single-cell datasets and large analysis results

---

## Quick Start

Because this is a multi-service research platform, each module may have its own dependencies and runtime assumptions. A general workflow is:

```bash
git clone https://github.com/ouev5/panhorm.git
cd panhorm
```

Then enter the module you want to run, for example:

```bash
cd animal_hormone_8000
cp .env.example .env
```

Fill in local environment variables in `.env` before running the service.

Typical variables may include:

```text
SECRET_KEY
DATABASE_NAME
DATABASE_USER
DATABASE_PASSWORD
DATABASE_HOST
DEEPSEEK_API_KEY
SERPER_API_KEY
OPENAI_API_KEY
```

> Never commit `.env` files or real credentials to GitHub.

---

## Data Restoration Notes

Some modules require large external data files to reproduce the full production behavior. These files are intentionally not included in the repository.

Examples:

- **AHormoneDB**: BLAST databases, RAG embeddings, text chunks, metadata, local embedding models, uploaded media files
- **GeneReg**: ENCODE ChIP-seq peak files under `encode_peaks/`
- **Venn / Single-cell tools**: single-cell datasets such as `.h5ad`, `.zarr`, or other large expression matrices
- **Multi-omics platform**: local data and generated result directories

Recommended storage options:

- Institutional storage
- Object storage
- Private data archive
- Separate private repository or release asset if appropriate

---

## Security Policy

This repository was sanitized before publication:

- Real API keys and service credentials were removed
- Runtime databases and user-uploaded files were removed
- Backup files and logs were removed
- Large private datasets were excluded
- Placeholder values such as `YOUR_API_KEY_HERE` are used where configuration is required

Before deploying or making substantial changes, run your own secret scan, for example:

```bash
find . -type f \( -name '.env' -o -name '*.key' -o -name '*.pem' -o -name '*.sqlite3' -o -name '*.db' \)
git grep -n -i -E "(SECRET_KEY|api[_-]?key|token|password|passwd|private key)" || true
```

For public releases, tools such as `gitleaks` or `trufflehog` are also recommended.

---

## Suggested Citation

If you use this codebase or refer to the platform in academic work, please cite this repository:

```text
ouev5. PanHorm Bioinformatics Platform: an integrated codebase for hormone-related bioinformatics analysis, gene regulation prediction, and multi-omics exploration. GitHub repository, 2026. https://github.com/ouev5/panhorm
```

A machine-readable citation file is also provided in [`CITATION.cff`](./CITATION.cff).

---

## License

This project is released under the **MIT License**. See the [`LICENSE`](./LICENSE) file for details.

---

## Maintainer

Maintained by **ouev5**.

Repository: <https://github.com/ouev5/panhorm>
