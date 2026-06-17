# Bioinformatics Analysis Platform: Technical Documentation

## Gene Insights Hub & AI Customer Service System

**Inner Mongolia University of Science and Technology**  
**School of Life Science and Technology**

---

## Abstract

This document describes the technical architecture and implementation of an integrated bioinformatics analysis platform consisting of two main components: (1) Gene Insights Hub, a multi-dimensional gene set analysis tool featuring Venn diagrams, UpSet plots, GO/KEGG enrichment analysis, single-cell expression visualization, and Sankey-bubble plots; and (2) an AI-powered customer service system utilizing Retrieval-Augmented Generation (RAG) with DeepSeek API for intelligent Q&A support. The platform is designed for researchers in molecular biology, genomics, and bioinformatics, providing accessible tools for gene set intersection analysis, functional enrichment, and knowledge retrieval.

---

## 1. Introduction

### 1.1 Background

Modern genomic studies frequently involve the analysis of multiple gene sets derived from differential expression analysis, genome-wide association studies, or literature mining. Visualizing intersections between gene sets and performing functional enrichment analysis are essential steps in biological interpretation. Traditional Venn diagrams become impractical when comparing more than 4-5 gene sets, and researchers often need to switch between multiple tools for different analysis steps.

### 1.2 Objectives

We developed an integrated platform to:
1. Provide intuitive visualization of multi-set gene intersections using UpSet plots
2. Enable seamless functional enrichment analysis through GO and KEGG databases
3. Offer single-cell gene expression visualization using real transcriptomic data
4. Support intelligent Q&A through an AI assistant trained on platform-specific knowledge

---

## 2. System Architecture

### 2.1 Overview

The platform adopts a microservices architecture with two independent services:

```
┌─────────────────────────────────────────────────────────────┐
│                    User Interface Layer                      │
│  ┌─────────────────┐  ┌─────────────────┐  ┌──────────────┐ │
│  │  Gene Insights  │  │  Animal Hormone │  │ AI Assistant │ │
│  │      Hub        │  │    Database     │  │   (Chat)     │ │
│  │  (Port 8001)    │  │   (Port 8000)   │  │ (Port 8002)  │ │
│  └─────────────────┘  └─────────────────┘  └──────────────┘ │
└─────────────────────────────────────────────────────────────┘
                              │
┌─────────────────────────────────────────────────────────────┐
│                    Application Layer                         │
│  ┌─────────────────────────────────────────────────────────┐│
│  │         Flask Application (Gene Insights Hub)           ││
│  │  - Venn/UpSet Analysis  - GO/KEGG Enrichment            ││
│  │  - Single-Cell Query    - Sankey-Bubble Plot            ││
│  └─────────────────────────────────────────────────────────┘│
│  ┌─────────────────────────────────────────────────────────┐│
│  │       FastAPI Application (AI Customer Service)          ││
│  │  - RAG Retrieval        - DeepSeek LLM Integration      ││
│  │  - Knowledge Base (812 Q&A pairs)                       ││
│  └─────────────────────────────────────────────────────────┘│
└─────────────────────────────────────────────────────────────┘
                              │
┌─────────────────────────────────────────────────────────────┐
│                      Data Layer                              │
│  ┌────────────┐  ┌────────────┐  ┌────────────────────────┐ │
│  │   MySQL    │  │  Enrichr   │  │   Tabula Sapiens       │ │
│  │  Database  │  │    API     │  │   Single-Cell Data     │ │
│  │            │  │            │  │  (94,836 cells)        │ │
│  └────────────┘  └────────────┘  └────────────────────────┘ │
│  ┌────────────┐  ┌────────────┐  ┌────────────────────────┐ │
│  │  STRING    │  │  Embedding │  │    Knowledge Base      │ │
│  │   API      │  │   Model    │  │    (812 entries)       │ │
│  │            │  │ (bge-large)│  │                        │ │
│  └────────────┘  └────────────┘  └────────────────────────┘ │
└─────────────────────────────────────────────────────────────┘
```

### 2.2 Technology Stack

| Component | Technology | Version |
|-----------|------------|---------|
| **Backend Framework** | Flask | 3.x |
| **API Framework** | FastAPI | Latest |
| **Database** | MySQL | 8.0 |
| **Python** | CPython | 3.12 |
| **Statistical Computing** | R | 4.3.3 |
| **Single-Cell Analysis** | Scanpy, scipy | Latest |
| **Visualization** | ggplot2, ggalluvial (R) | Latest |
| **Embedding Model** | BAAI/bge-large-zh-v1.5 | Latest |
| **LLM API** | DeepSeek Chat | deepseek-chat |
| **Deployment** | systemd, Gunicorn | Latest |

---

## 3. Gene Insights Hub Implementation

### 3.1 Core Features

#### 3.1.1 Multi-Set Venn/UpSet Analysis

**Implementation Details:**
- Users input gene lists through text boxes, file upload, or database search
- Intersection calculation uses Python's `itertools.combinations` for efficient set operations
- Results displayed as UpSet plots (matrix view) and classic Venn diagrams
- Interactive selection: clicking bars selects corresponding gene subsets

**Technical Approach:**
```python
from itertools import combinations

def calculate_intersections(gene_lists):
    """Calculate all intersections between gene lists"""
    list_names = list(gene_lists.keys())
    results = {}
    
    for r in range(1, len(list_names) + 1):
        for combo in combinations(list_names, r):
            intersection = set(gene_lists[combo[0]])
            for name in combo[1:]:
                intersection &= set(gene_lists[name])
            
            # Find genes unique to this intersection
            exclusive = intersection.copy()
            for name in list_names:
                if name not in combo:
                    exclusive -= set(gene_lists[name])
            
            results[combo] = exclusive
    
    return results
```

#### 3.1.2 GO/KEGG Enrichment Analysis

**Data Source:** Enrichr API (https://maayanlab.cloud/Enrichr/)

**Gene Set Libraries:**
- GO_Biological_Process_2023: Gene Ontology biological processes
- KEGG_2021_Human: KEGG metabolic and signaling pathways

**Implementation:**
```python
import requests

def enrichment_analysis(genes, analysis_type="go"):
    """Perform enrichment analysis via Enrichr API"""
    # Select gene set library
    gene_library = "GO_Biological_Process_2023" if analysis_type == "go" else "KEGG_2021_Human"
    
    # Enrichr API endpoints
    ENRICHR_URL = "https://maayanlab.cloud/Enrichr/addList"
    QUERY_URL = "https://maayanlab.cloud/Enrichr/enrich"
    
    # Submit gene list
    genes_str = "\n".join(genes)
    response = requests.post(ENRICHR_URL, files={"list": (None, genes_str)})
    user_list_id = response.json()["userListId"]
    
    # Query enrichment results
    query_response = requests.get(f"{QUERY_URL}?userListId={user_list_id}&backgroundType={gene_library}")
    results = query_response.json()[gene_library]
    
    # Parse results: term, p-value, adjusted p-value, genes, combined score
    return [{
        "term": r[1],
        "pvalue": r[2],
        "adj_pvalue": r[6],
        "genes": r[5],
        "combined_score": r[4]
    } for r in results[:25]]
```

**Output Formats:**
- Bar chart: Top enriched terms by combined score
- Bubble plot: Enrichment ratio vs. significance
- Pie chart: Proportional distribution
- Data table: Complete statistics for download

#### 3.1.3 Single-Cell Expression Visualization

**Data Source:** Tabula Sapiens (Human Cell Atlas)

**Dataset Characteristics:**
- **Cells:** 94,836 single cells
- **Genes:** 34,562 genes
- **Tissues:** 24 tissue types
- **Cell Types:** 161 cell types

**Reference:** Tabula Sapiens Consortium. "The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans." *Science* 376, no. 6584 (2022): eabl4896.

**Implementation:**
```python
import scanpy as sc
import numpy as np
from scipy import sparse

class SingleCellQuery:
    def __init__(self, data_path):
        # Load pre-processed data
        self.expr_matrix = sparse.load_npz(f"{data_path}/expr_matrix.npz")
        self.gene_list = json.load(open(f"{data_path}/gene_list.json"))
        self.cell_metadata = json.load(open(f"{data_path}/cell_metadata.json"))
        
        # Build gene index
        self.gene_to_idx = {g: i for i, g in enumerate(self.gene_list)}
    
    def query_gene_expression(self, gene_name):
        """Query expression of a gene across all cells"""
        if gene_name not in self.gene_to_idx:
            return None
        
        gene_idx = self.gene_to_idx[gene_name]
        expr = self.expr_matrix[gene_idx, :].toarray().flatten()
        
        # Return expression values with cell metadata
        return {
            "expression": expr.tolist(),
            "umap_coords": self.umap_coords,
            "cell_types": self.cell_types,
            "tissues": self.tissues
        }
```

**Visualization Methods:**
- **UMAP:** Pre-computed coordinates for fast rendering
- **t-SNE:** Real-time computation using sklearn.manifold.TSNE
- **Color encoding:** Cell type, tissue origin, or gene expression level

#### 3.1.4 Sankey-Bubble Enrichment Plot

**Purpose:** Visualize gene-pathway relationships in a combined view

**Implementation (R/ggplot2):**
```r
library(ggplot2)
library(ggalluvial)

generate_sankey_bubble <- function(enrichment_data, output_file) {
  # Sankey plot: gene-to-pathway flow
  sankey <- ggplot(data, aes(x = axis, stratum = stratum, 
                              alluvium = alluvium, fill = fill)) +
    geom_alluvium(width = 0.3) +
    geom_stratum(width = 0.3) +
    theme_minimal()
  
  # Bubble plot: enrichment statistics
  bubble <- ggplot(data, aes(x = ratio, y = term, 
                              size = gene_count, color = -log10(pvalue))) +
    geom_point() +
    scale_size_continuous(range = c(3, 10)) +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal()
  
  # Combine plots
  combined <- sankey + bubble + plot_layout(ncol = 2)
  ggsave(output_file, combined, width = 14, height = 10, dpi = 300)
}
```

**Output:** Publication-ready figures (PNG/JPG/SVG)

#### 3.1.5 STRING Network Integration

**Purpose:** Protein-protein interaction network visualization

**Implementation:**
```python
def generate_string_url(genes, species=9606):
    """Generate STRING database URL for protein network"""
    base_url = "https://string-db.org/cgi/network"
    gene_str = "%0d".join(genes[:50])  # Limit to 50 genes
    
    return f"{base_url}?identifiers={gene_str}&species={species}"
```

**Features:**
- Directs users to STRING database for network visualization
- Supports up to 50 genes per query
- Species: Human (9606)

### 3.2 Data Preprocessing Toolkit

**Available Functions:**
1. **Deduplication:** Remove duplicate gene names
2. **Format Standardization:** Convert to official gene symbols
3. **Case Conversion:** Uppercase/lowercase normalization
4. **Empty Line Removal:** Clean whitespace
5. **Annotation Removal:** Strip gene descriptions

### 3.3 Deployment Configuration

**Systemd Service Unit:**
```ini
[Unit]
Description=Venn Tool - Gene Analysis Service
After=network.target

[Service]
Type=simple
User=root
WorkingDirectory=/www/wwwroot/venn-tool
Environment="PATH=/www/wwwroot/venv/bin"
ExecStart=/www/wwwroot/venv/bin/python3 app.py
Restart=always
RestartSec=5

[Install]
WantedBy=multi-user.target
```

**Port:** 8001  
**URL:** http://43.99.62.219:8001/

---

## 4. AI Customer Service System

### 4.1 Architecture

The AI customer service system employs a **Retrieval-Augmented Generation (RAG)** architecture combining semantic search with large language model generation.

```
┌─────────────────────────────────────────────────────────────┐
│                      User Query                              │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│              Knowledge Base Retrieval                        │
│  ┌─────────────────────────────────────────────────────────┐│
│  │  1. Query Embedding (BAAI/bge-large-zh-v1.5)            ││
│  │  2. Vector Similarity Search (Cosine Similarity)        ││
│  │  3. Top-K Relevant Passages Retrieval                   ││
│  └─────────────────────────────────────────────────────────┘│
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│                  Context Construction                        │
│  System Prompt + Retrieved Knowledge + User Query           │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│              DeepSeek LLM Generation                         │
│  Model: deepseek-chat                                       │
│  API: https://api.deepseek.com/v1/chat/completions          │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│                     Response Output                          │
│  Markdown-formatted answer with code blocks, lists, etc.    │
└─────────────────────────────────────────────────────────────┘
```

### 4.2 Knowledge Base

**Content Categories:**
| Category | Description | Entries |
|----------|-------------|---------|
| Platform Functions | Feature descriptions and usage | 150+ |
| GO/KEGG Analysis | Enrichment analysis methods | 80+ |
| Single-Cell Analysis | UMAP/t-SNE visualization | 50+ |
| Hormone Knowledge | Biological mechanisms | 100+ |
| Troubleshooting | Common issues and solutions | 60+ |
| Operation Guides | Step-by-step instructions | 200+ |
| Database Knowledge | BLAST, STRING, KEGG | 70+ |
| Workflow Guides | Complete analysis pipelines | 50+ |

**Total Entries:** 812 Q&A pairs + 10 knowledge documents

### 4.3 Embedding Model

**Model:** BAAI/bge-large-zh-v1.5

**Specifications:**
- Dimension: 1024
- Max Sequence Length: 512
- Language: Chinese & English bilingual
- Performance: C-MTEB benchmark leader

**Implementation:**
```python
from sentence_transformers import SentenceTransformer
import numpy as np

class KnowledgeRetriever:
    def __init__(self, model_name="BAAI/bge-large-zh-v1.5"):
        self.model = SentenceTransformer(model_name)
        self.knowledge_base = self.load_knowledge()
        self.embeddings = self.encode_knowledge()
    
    def encode_knowledge(self):
        """Encode all knowledge entries"""
        texts = [f"{k['instruction']} {k['output']}" 
                 for k in self.knowledge_base]
        return self.model.encode(texts, normalize_embeddings=True)
    
    def retrieve(self, query, top_k=5):
        """Retrieve relevant knowledge for query"""
        query_embedding = self.model.encode([query], normalize_embeddings=True)[0]
        
        # Cosine similarity
        similarities = np.dot(self.embeddings, query_embedding)
        top_indices = np.argsort(similarities)[-top_k:][::-1]
        
        return [self.knowledge_base[i] for i in top_indices]
```

### 4.4 DeepSeek Integration

**API Configuration:**
```python
DEEPSEEK_API_URL = "https://api.deepseek.com/v1/chat/completions"
MODEL = "deepseek-chat"

def generate_response(query, context):
    """Generate response using DeepSeek API"""
    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": f"Context: {context}\n\nQuestion: {query}"}
    ]
    
    response = requests.post(
        DEEPSEEK_API_URL,
        headers={
            "Authorization": f"Bearer {DEEPSEEK_API_KEY}",
            "Content-Type": "application/json"
        },
        json={
            "model": MODEL,
            "messages": messages,
            "temperature": 0.7,
            "max_tokens": 2000
        }
    )
    
    return response.json()["choices"][0]["message"]["content"]
```

### 4.5 System Prompt

```
You are an AI assistant for a bioinformatics analysis platform developed by 
Inner Mongolia University of Science and Technology.

The platform contains two main modules:
1. Gene Insights Hub (Venn Diagram Analysis Tool) - Supports UpSet plots, 
   GO enrichment analysis, KEGG analysis, single-cell analysis, sankey bubble 
   plots, and STRING network visualization.
2. Animal Hormone Regulation Database - Hormone data query, BLAST alignment, 
   transcriptome analysis, DNA methylation analysis, GNN prediction, 
   phylogenetic analysis, and KEGG pathway analysis.

Your responsibilities:
- Answer questions about platform functions enthusiastically and friendly
- Explain bioinformatics analysis methods
- Provide hormone-related knowledge
- Guide users through operation steps
- Answer common questions

Response requirements:
- Professional, accurate, and friendly
- Answer in ENGLISH only
- Use clear structure with bullet points when appropriate
- Provide detailed steps for specific operations
```

### 4.6 Frontend Integration

**Integration Points:**
1. **Floating Button:** Right-bottom corner on all pages
2. **Chat Window:** Pop-up dialog with quick questions
3. **Markdown Rendering:** Support for formatted responses

**Quick Question Buttons:**
- Platform Introduction
- GO Analysis
- BLAST
- Single Cell Analysis

---

## 5. Database Integration

### 5.1 MySQL Database

**Connection Configuration:**
```python
import pymysql

def get_database_connection():
    return pymysql.connect(
        host='127.0.0.1',
        user='hhh',
        database='animal_hormone',
        charset='utf8mb4',
        cursorclass=pymysql.cursors.DictCursor
    )
```

**Tables Used:**
- `hormone_related_genes`: Gene-hormone associations
- `hormone_receptor_full`: Receptor information
- `hormone_receptor_info`: Additional receptor metadata

### 5.2 External APIs

| API | Purpose | Rate Limit |
|-----|---------|------------|
| Enrichr | GO/KEGG enrichment | ~1 request/sec |
| STRING | Protein interactions | No limit |
| DeepSeek | LLM generation | API quota |

---

## 6. User Interface Design

### 6.1 Design Principles

1. **Responsive Design:** Mobile-first approach with Tailwind CSS
2. **Progressive Disclosure:** Complex features in expandable sections
3. **Visual Feedback:** Loading indicators and success/error messages
4. **Accessibility:** High contrast, clear typography, keyboard navigation

### 6.2 Color Scheme

| Element | Color | Hex Code |
|---------|-------|----------|
| Primary | Blue | #165DFF |
| Secondary | Teal | #36CFC9 |
| Accent | Purple | #722ED1 |
| Success | Green | #52C41A |
| Warning | Amber | #FAAD14 |
| Danger | Red | #FF4D4F |

### 6.3 Interactive Features

- **Click-to-select:** UpSet plot bars for gene selection
- **Real-time filtering:** Database search with autocomplete
- **Drag-and-drop:** File upload for gene lists
- **Hover tooltips:** Additional information on demand

---

## 7. Performance Optimization

### 7.1 Single-Cell Data

**Challenge:** 94,836 cells × 34,562 genes = 3.28 billion data points

**Solutions:**
1. **Sparse Matrix Storage:** scipy.sparse.csr_matrix (98% space reduction)
2. **Pre-computed UMAP:** Avoid real-time dimensionality reduction
3. **Lazy Loading:** Load data on-demand
4. **Caching:** Cache frequently queried genes

**Performance Metrics:**
- Data loading: ~2 seconds
- Single gene query: <100ms
- t-SNE computation: ~30 seconds

### 7.2 API Optimization

1. **Connection Pooling:** Reuse database connections
2. **Request Batching:** Combine multiple enrichment queries
3. **Async Processing:** Non-blocking I/O for file operations
4. **Response Compression:** gzip for large JSON responses

---

## 8. Deployment and Operations

### 8.1 Server Specifications

| Resource | Specification |
|----------|---------------|
| CPU | Multi-core |
| RAM | Sufficient for single-cell data (~12GB for matrix) |
| Storage | SSD recommended |
| OS | Linux (Ubuntu/CentOS) |

### 8.2 Service Management

```bash
# Start services
systemctl start venn-tool
systemctl start bioinfo-customer-service

# Check status
systemctl status venn-tool
systemctl status bioinfo-customer-service

# View logs
journalctl -u venn-tool -f
journalctl -u bioinfo-customer-service -f

# Restart services
systemctl restart venn-tool
systemctl restart bioinfo-customer-service
```

### 8.3 Monitoring

- **Service Health:** HTTP endpoints for health checks
- **Log Aggregation:** journalctl for centralized logging
- **Error Tracking:** Structured logging with severity levels

---

## 9. Security Considerations

### 9.1 API Security

- API keys stored as environment variables
- Input validation on all endpoints
- Rate limiting for external API calls
- CORS configuration for cross-origin requests

### 9.2 Data Security

- Database credentials not exposed to frontend
- User data not stored on server
- HTTPS for all communications
- Regular security updates

---

## 10. Future Development

### 10.1 Planned Features

1. **User Accounts:** Save analysis history and gene lists
2. **Batch Analysis:** Process multiple gene sets simultaneously
3. **Custom Gene Sets:** User-defined gene set libraries
4. **Export Integration:** Direct export to publication formats
5. **Mobile App:** Native mobile application

### 10.2 Scalability

- Containerization with Docker
- Load balancing for multiple instances
- Database sharding for large datasets
- CDN for static assets

---

## 11. Conclusion

This integrated bioinformatics platform provides researchers with comprehensive tools for gene set analysis, from intersection visualization to functional enrichment and single-cell expression exploration. The AI customer service system enhances user experience by providing intelligent, context-aware assistance. The modular architecture allows for easy extension and maintenance, while the choice of open-source technologies ensures reproducibility and accessibility.

---

## References

1. Tabula Sapiens Consortium. "The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans." *Science* 376, no. 6584 (2022): eabl4896.

2. Chen, E.Y., et al. "Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool." *BMC Bioinformatics* 14, no. 1 (2013): 128.

3. Szklarczyk, D., et al. "STRING v11: protein-protein association networks with increased coverage, supporting functional discovery in genome-wide experimental datasets." *Nucleic Acids Research* 47, no. D1 (2019): D607-D613.

4. Xiao, S., et al. "C-Pack: Packaged C/C++ library for embedding-based search." *arXiv preprint arXiv:2402.16942* (2024).

5. Lex, A., et al. "UpSet: visualization of intersecting sets." *IEEE Transactions on Visualization and Computer Graphics* 20, no. 12 (2014): 1983-1992.

---

## Appendix A: API Endpoints

### Gene Insights Hub (Port 8001)

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/` | GET | Main page |
| `/api/parse_file` | POST | Parse uploaded gene file |
| `/api/analyze` | POST | Calculate gene intersections |
| `/api/enrichment` | POST | GO/KEGG enrichment analysis |
| `/api/sankey_bubble_plot` | POST | Generate Sankey-bubble plot |
| `/api/single_cell_query` | POST | Single-cell expression query |
| `/api/single_cell_query_tsne` | POST | t-SNE visualization |
| `/api/string_url` | POST | Generate STRING URL |
| `/api/search_database` | GET | Search gene database |

### AI Customer Service (Port 8002)

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/` | GET | Chat interface |
| `/api/chat` | POST | Chat with AI assistant |
| `/api/health` | GET | Service health check |
| `/docs` | GET | API documentation |

---

## Appendix B: File Structure

```
/www/wwwroot/venn-tool/
├── app.py                 # Main Flask application
├── single_cell_real.py    # Single-cell query module
├── tsne_api.py           # t-SNE API module
├── templates/
│   ├── index.html        # Main page template
│   └── ...
├── static/
│   ├── css/
│   └── js/
├── r_scripts/
│   └── sankey_bubble.R   # R script for Sankey plots
├── single_cell_data/
│   ├── expr_matrix.npz   # Sparse expression matrix
│   ├── gene_list.json    # Gene names
│   └── cell_metadata.json # Cell annotations
└── fine_tuning_data/
    ├── basic_qa.json     # Basic Q&A pairs
    ├── hormone_qa.json   # Hormone knowledge
    └── ...               # Other knowledge files

/www/wwwroot/customer_service/
├── bioinfo_customer_service.py  # FastAPI application
└── ...
```

---

## Appendix C: Code Availability

The platform is deployed at:
- Gene Insights Hub: http://43.99.62.219:8001/
- Animal Hormone Database: https://43.99.62.219/
- AI Assistant: http://43.99.62.219:8002/

Source code available upon request from the authors.

---

**Document Version:** 1.0  
**Last Updated:** April 19, 2026  
**Authors:** Inner Mongolia University of Science and Technology, School of Life Science and Technology
