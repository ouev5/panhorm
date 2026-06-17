from django.shortcuts import render, redirect, get_object_or_404
from django.http import JsonResponse, HttpResponse
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.models import User
from django.contrib.auth.decorators import login_required, user_passes_test
from django.utils import timezone
from datetime import datetime, timedelta
from django.core.paginator import Paginator, PageNotAnInteger, EmptyPage
from django.core.cache import cache
from django.db.models import Q, Count
from django.views.decorators.csrf import csrf_protect, csrf_exempt
from django.conf import settings
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from django.views import View
from django.utils.decorators import method_decorator
from django.apps import apps
from django.db.models import Model
from django.core.files.base import ContentFile
import requests
import re
import logging
import json
import os
import uuid
import tempfile
import zipfile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from weasyprint import HTML
import traceback
from urllib.parse import quote
import traceback
from xhtml2pdf import pisa
from django.template.loader import get_template
from io import BytesIO, StringIO
# RAG：

from .utils.rag_retriever import hybrid_advanced_retrieval, build_biogpt_prompt, vector_store, initialize_vector_store_async
from .utils import rag_deepseek 
from .utils.rag_deepseek import extract_entities_for_pdf
import json
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt


def network_data(request):
    """API: Return hormone-gene-disease network data for ECharts graph visualization"""
    from django.http import JsonResponse
    from hormone_app.models import HormoneRelatedGene
    from collections import defaultdict
    
    query = request.GET.get('q', '').strip()
    limit = int(request.GET.get('limit', '5'))  # Number of hormones to show (default 5)
    max_hormone_records = int(request.GET.get('max_records', '15'))  # Max records per hormone
    
    # If search query provided, search by hormone name
    if query:
        hormones = HormoneRelatedGene.objects.filter(
            Q(hormone_name__icontains=query) | Q(related_genes__icontains=query) | Q(related_diseases__icontains=query)
        ).values('hormone_name').annotate(
            count=Count('id')
        ).order_by('-count')[:limit]
    else:
        # Default: top hormones by record count
        hormones = HormoneRelatedGene.objects.values('hormone_name').annotate(
            count=Count('id')
        ).order_by('-count')[:limit]
    
    nodes = {}  # name -> node data
    links = []
    
    for h in hormones:
        hormone_name = h['hormone_name']
        records = HormoneRelatedGene.objects.filter(
            hormone_name=hormone_name
        ).exclude(
            related_genes='Not available'
        ).exclude(
            related_genes=''
        ).values('related_genes', 'related_diseases', 'regulation_type', 'organism')[:max_hormone_records]
        
        # Add hormone node
        if hormone_name not in nodes:
            nodes[hormone_name] = {
                'name': hormone_name,
                'category': 0,
                'symbolSize': 65,
            }
        
        for rec in records:
            gene = rec.get('related_genes', '')
            disease = rec.get('related_diseases', '')
            reg_type = rec.get('regulation_type', '')
            organism = rec.get('organism', '')
            
            # Handle multiple genes in one record (comma-separated)
            genes = [g.strip() for g in gene.split(',') if g.strip()]
            
            for g in genes:
                g_name = g.split(' (')[0].strip()  # Remove parenthetical descriptions
                if not g_name or len(g_name) > 30:
                    continue
                
                if g_name not in nodes:
                    nodes[g_name] = {
                        'name': g_name,
                        'category': 1,
                        'symbolSize': 45,
                    }
                
                # Hormone -> Gene link
                link_label = reg_type if reg_type else ''
                links.append({
                    'source': hormone_name,
                    'target': g_name,
                    'value': link_label,
                })
                
                # Gene -> Disease link
                if disease and disease != 'Not available' and len(disease) < 50:
                    diseases = [d.strip() for d in disease.split(';') if d.strip()]
                    for dis in diseases[:2]:  # Max 2 diseases per gene
                        dis_name = dis[:40]  # Truncate long names
                        if dis_name not in nodes:
                            nodes[dis_name] = {
                                'name': dis_name,
                                'category': 2,
                                'symbolSize': 35,
                            }
                        links.append({
                            'source': g_name,
                            'target': dis_name,
                            'value': '',
                        })
    
    return JsonResponse({
        'nodes': list(nodes.values()),
        'links': links,
        'categories': [
            {'name': 'Hormone'},
            {'name': 'Gene'},
            {'name': 'Disease'},
        ],
    })


from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer  # 按需导入其他组件
from reportlab.lib.styles import getSampleStyleSheet  # 可能用到的样式库
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image  # 补充 Image
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib import colors
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import inch
from urllib.parse import quote  # 用于处理文件名编码
import logging
# 配置日志
logger = logging.getLogger(__name__)

# BioGPT API地址（同一服务器）
BIOGPT_API_URL = "http://localhost:8001/generate"
# 预定义基础术语（与biogpt_api.py保持一致）
BASIC_DEFINITIONS = {
    "hormone": "A hormone is a chemical messenger produced by endocrine glands that travels through the bloodstream to regulate the activity of target cells or organs. Hormones control essential bodily functions including growth, metabolism, reproduction, and homeostasis.",
    "gene": "A gene is a segment of DNA that contains the instructions for building one or more molecules that help the body work. Genes are the basic units of heredity and determine traits in living organisms.",
    "protein": "A protein is a large biomolecule made of amino acids that performs a vast array of functions in living organisms, including catalyzing metabolic reactions, DNA replication, and providing structural support to cells."
}

# 导入模型和表单
from .models import (
    Hormone, Literature, AnalysisHistory, HormoneRawData,
    HormoneRelatedGene, HormoneReceptorInfo, HormoneReceptorFull,
    GeneAnalysisTask, GeneUploadedFile, GeneAnalysisResult, AnalysisTask, UploadedFile, 
    AnalysisParameter, AnalysisResult, AnalysisJob
)
from .forms import (
    AnalysisJobForm, FileUploadForm, QualityControlForm,
    PreprocessingForm, DifferentialMethylationForm,
    RegionAnalysisForm, CorrelationForm, EnrichmentForm
)
from .tasks import run_analysis_pipeline  # 异步任务
from .biogpt_utils import generate_biomed_response
import urllib.parse
# ===================== PDF报告生成视图 =====================
# 辅助函数：将HTML转换为PDF
def render_to_pdf(template_src, context_dict={}):
    template = get_template(template_src)
    html = template.render(context_dict)
    result = BytesIO()
    # 生成PDF（指定中文字体支持，避免乱码）
    pdf = pisa.pisaDocument(
        BytesIO(html.encode("UTF-8")),
        result,
        encoding="UTF-8",
        # 若需要支持中文，需指定中文字体路径（可选）
        # default_css="""
        #     @font-face {
        #         font-family: SimHei;
        #         src: url("/static/fonts/simhei.ttf");
        #     }
        #     body { font-family: SimHei, Arial; }
        # """
    )
    if not pdf.err:
        return HttpResponse(result.getvalue(), content_type='application/pdf')
    return None

from django.shortcuts import render, redirect, get_object_or_404
from django.http import JsonResponse, HttpResponse
from django.utils import timezone
from django.db.models import Q, Count
import requests
import re
import logging
import urllib.parse
# 确保导入相关模型（HormoneRelatedGene、HormoneReceptorFull、HormoneReceptorInfo）

logger = logging.getLogger(__name__)
BIOGPT_API_URL = "http://localhost:8001/generate"

def generate_hormone_report(request):
    """生成激素/基因/疾病相关的PDF报告并提供下载（支持精确/模糊匹配，完全兼容大小写/格式差异）"""
    # 1. 获取查询参数（新增 disease 参数，支持三种查询场景）
    hormone_name = request.GET.get("hormone", "").strip()
    raw_gene_query = request.GET.get("gene", "").strip()
    disease_name = request.GET.get("disease", "").strip()  # 新增：疾病查询参数
    match_type = request.GET.get("match_type", "fuzzy").lower()
    
    # 校验参数
    if match_type not in ["exact", "fuzzy"]:
        return JsonResponse({
            "code": 400,
            "msg": "match_type参数值无效，仅支持exact（精确匹配）或fuzzy（模糊匹配）"
        })
    # 新增：校验三种查询参数至少传入一个
    if not hormone_name and not raw_gene_query and not disease_name:
        return JsonResponse({
            "code": 400,
            "msg": "请传入激素名称（?hormone=xxx）、基因名称（?gene=xxx）或疾病名称（?disease=xxx）"
        })
    
    # 2. 初始化数据容器
    related_genes = []
    receptor_full = []
    receptor_info = []
    total_records = 0
    summary = ""

    # 优化：文本清洗函数（保留有效数据格式）
    def clean_text(text):
        if not text or str(text).strip() in ["", "N/A"]:
            return "N/A"
        text_str = str(text).strip()
        text_str = re.sub(r'<script.*?</script>', ' ', text_str, flags=re.DOTALL)
        text_str = re.sub(r'<[^>]+>', ' ', text_str)
        text_str = re.sub(r'[\r\n]{2,}', '\n', text_str)
        text_str = re.sub(r'\s{2,}', ' ', text_str)
        return text_str
    
    # 核心：基因名称标准化（兼容大小写、分隔符）
    def standardize_gene_name(gene_name):
        if not gene_name:
            return ""
        standardized = gene_name.strip().upper()  # 转为大写
        standardized = re.sub(r'[-_.]', '', standardized)  # 去除分隔符
        standardized = re.sub(r'\s+', '', standardized)  # 去除空格
        return standardized
    
    # 优化：激素名称标准化（兼容大小写）
    def standardize_hormone_name(hormone_name):
        if not hormone_name:
            return ""
        # 首字母大写，其余小写（兼容数据库中首字母大写的格式，如Glucagon）
        return hormone_name.strip().capitalize()
    
    # 新增：疾病名称标准化（兼容大小写，去除冗余空格）
    def standardize_disease_name(disease_name):
        if not disease_name:
            return ""
        # 首字母大写，其余小写，去除多余空格（如 "type 2 diabetes" → "Type 2 diabetes"）
        disease_str = disease_name.strip()
        if not disease_str:
            return ""
        # 拆分单词，首字母大写（保留数字和特殊格式）
        words = disease_str.split()
        standardized_words = []
        for word in words:
            if word.isdigit() or (len(word) > 0 and word[0].isdigit()):
                standardized_words.append(word)
            else:
                standardized_words.append(word.capitalize())
        return " ".join(standardized_words)
    
    # 核心：查询函数（完全兼容大小写/格式差异，新增疾病字段支持）
    def query_records(model, field, value, match_type):
        """
        通用查询函数：
        - 基因字段：标准化后匹配（大写+去分隔符）
        - 激素字段：大小写不敏感匹配
        - 疾病字段：大小写不敏感匹配（新增）
        - 其他字段：保留原有逻辑
        """
        gene_fields = ["related_genes", "receptor_coding_genes"]
        hormone_fields = ["hormone_name"]
        disease_fields = ["related_diseases"]  # 新增：疾病相关字段（对应HormoneRelatedGene的related_diseases）
        is_gene_field = field in gene_fields
        is_hormone_field = field in hormone_fields
        is_disease_field = field in disease_fields
        
        # 标准化查询值
        query_value = ""
        if is_gene_field:
            query_value = standardize_gene_name(value)
            backup_query_value = value.strip().upper()  # 备份原始大写值
        elif is_hormone_field:
            query_value = value.strip()  # 激素字段不改变原始值，用icontains兼容大小写
        elif is_disease_field:
            query_value = value.strip()  # 疾病字段不改变原始值，用icontains兼容大小写
        else:
            query_value = value.strip()
        
        if not query_value:
            return model.objects.none()
        
        # 执行查询
        if match_type == "exact":
            if is_gene_field:
                # 基因精确匹配：同时匹配标准化值和原始大写值
                return model.objects.filter(
                    Q(**{f"{field}__exact": query_value}) |
                    Q(**{f"{field}__exact": backup_query_value})
                ).all()
            elif is_hormone_field or is_disease_field:
                # 激素/疾病精确匹配：大小写不敏感（兼容数据库大小写）
                return model.objects.filter(**{f"{field}__iexact": query_value}).all()
            else:
                return model.objects.filter(**{f"{field}__exact": query_value}).all()
        else:  # fuzzy
            if is_gene_field:
                # 基因模糊匹配：标准化值+原始大写值，双重匹配
                return model.objects.filter(
                    Q(**{f"{field}__icontains": query_value}) |
                    Q(**{f"{field}__icontains": backup_query_value})
                ).all()
            elif is_hormone_field or is_disease_field:
                # 激素/疾病模糊匹配：大小写不敏感（新增疾病字段适配）
                return model.objects.filter(**{f"{field}__icontains": query_value}).all()
            else:
                return model.objects.filter(**{f"{field}__icontains": query_value}).all()
    
    try:
        # 3. 查询相关基因数据（核心：支持激素/基因/疾病三种查询，兼容大小写）
        hg_records = []
        if hormone_name:
            # 激素查询（原有逻辑保留）
            hg_records = query_records(HormoneRelatedGene, "hormone_name", hormone_name, match_type)
        elif raw_gene_query:
            # 基因查询（原有逻辑保留）
            hg_records = query_records(HormoneRelatedGene, "related_genes", raw_gene_query, match_type)
        elif disease_name:
            # 新增：疾病查询（基于related_diseases字段）
            hg_records = query_records(HormoneRelatedGene, "related_diseases", disease_name, match_type)
        
        # 打印查询日志（调试用，可删除）
        logger.debug(f"查询类型：{'激素' if hormone_name else '基因' if raw_gene_query else '疾病'}")
        logger.debug(f"查询值：{hormone_name or raw_gene_query or disease_name}")
        logger.debug(f"标准化后（疾病）：{standardize_disease_name(disease_name)}")
        logger.debug(f"匹配到的相关基因记录数：{hg_records.count()}")
        
        # 整理相关基因数据
        for hg in hg_records:
            related_genes.append({
                "Hormone_Name": clean_text(hg.hormone_name),
                "Related_Genes": clean_text(hg.related_genes),
                "PMID": clean_text(hg.pmid),
                "Related_Diseases": clean_text(hg.related_diseases),
                "Gene_Sequence": "序列数据已隐藏" if hg.gene_sequence else "N/A"
            })
        
        # 4. 查询激素-受体完整信息（优化：支持疾病查询衍生的激素名称匹配）
        hrf_records = []
        # 提取所有有效激素名称（用于关联受体数据，无论哪种查询类型）
        valid_hormone_names = list({
            item["Hormone_Name"] for item in related_genes 
            if item["Hormone_Name"] != "N/A" and item["Hormone_Name"].strip()
        })
        
        if hormone_name:
            # 激素查询：直接匹配受体数据
            hrf_records = query_records(HormoneReceptorFull, "hormone_name", hormone_name, match_type)
        elif raw_gene_query or disease_name:
            # 基因/疾病查询：基于提取的有效激素名称，匹配受体数据
            logger.debug(f"从查询结果中提取的有效激素名称：{valid_hormone_names}")
            if valid_hormone_names:
                if match_type == "exact":
                    # 精确匹配：激素名称大小写不敏感
                    query = Q()
                    for name in valid_hormone_names:
                        query |= Q(hormone_name__iexact=name)
                    hrf_records = HormoneReceptorFull.objects.filter(query).all()
                else:
                    # 模糊匹配：激素名称大小写不敏感
                    query = Q()
                    for name in valid_hormone_names:
                        query |= Q(hormone_name__icontains=name)
                    hrf_records = HormoneReceptorFull.objects.filter(query).all()
        
        # 整理受体完整信息
        for hrf in hrf_records:
            receptor_full.append({
                "Hormone_Name": clean_text(hrf.hormone_name),
                "Hormone_Species_Name": clean_text(hrf.hormone_species_name),
                "Hormone_UniProt_ID": clean_text(hrf.hormone_uniprot_id),
                "Receptor_Name": clean_text(hrf.receptor_name),
                "Receptor_Species_Name": clean_text(hrf.receptor_species_name),
                "Receptor_Coding_Genes": clean_text(hrf.receptor_coding_genes)
            })
        
        # 5. 查询激素受体详细信息（与受体完整信息逻辑一致，支持三种查询类型）
        hri_records = []
        if hormone_name:
            hri_records = query_records(HormoneReceptorInfo, "hormone_name", hormone_name, match_type)
        elif raw_gene_query or disease_name:
            if valid_hormone_names:
                if match_type == "exact":
                    query = Q()
                    for name in valid_hormone_names:
                        query |= Q(hormone_name__iexact=name)
                    hri_records = HormoneReceptorInfo.objects.filter(query).all()
                else:
                    query = Q()
                    for name in valid_hormone_names:
                        query |= Q(hormone_name__icontains=name)
                    hri_records = HormoneReceptorInfo.objects.filter(query).all()
        
        # 整理受体详细信息
        for hri in hri_records:
            receptor_info.append({
                "Hormone_Name": clean_text(hri.hormone_name),
                "PubChem_ID": clean_text(hri.pubchem_id),
                "Receptor_Name": clean_text(hri.receptor_name),
                "Receptor_UniProt_ID": clean_text(hri.receptor_uniprot_id),
                "Receptor_Species_Name": clean_text(hri.receptor_species_name),
                "Receptor_Coding_Genes": clean_text(hri.receptor_coding_genes)
            })
        
        # 计算有效记录数
        real_related_genes = [g for g in related_genes if g["Hormone_Name"] != "N/A"]
        real_receptor_full = [r for r in receptor_full if r["Hormone_Name"] != "N/A"]
        real_receptor_info = [i for i in receptor_info if i["Hormone_Name"] != "N/A"]
        total_records = len(real_related_genes) + len(real_receptor_full) + len(real_receptor_info)
        logger.debug(f"有效记录总数：{total_records}")
        
        # 6. 调用BiogPT生成总结（基于有效数据，适配三种查询类型）
        if total_records > 0:
            # 确定查询目标
            query_type = "激素" if hormone_name else "基因" if raw_gene_query else "疾病"
            query_target = hormone_name or raw_gene_query or disease_name
            # 构造更精准的总结查询文本
            disease_list = [g["Related_Diseases"] for g in real_related_genes if g["Related_Diseases"] != "N/A"]
            hormone_list = [g["Hormone_Name"] for g in real_related_genes if g["Hormone_Name"] != "N/A"]
            query_text = f"分析{query_type}{query_target}的生物医学相关数据：相关激素为{list(set(hormone_list))}，相关疾病为{list(set(disease_list))}，生成专业的生物医学总结（1-2句话，简洁准确）"
            
            try:
                bio_response = requests.post(
                    BIOGPT_API_URL,
                    json={"input_text": query_text},
                    timeout=30
                )
                bio_response.raise_for_status()
                raw_summary = bio_response.json().get("response", "未获取到分析总结")
                summary = clean_text(raw_summary)
            except Exception as e:
                summary = f"BiogPT总结获取失败：{str(e)[:50]}"
                logger.warning(f"BiogPT调用失败：{e}")
        else:
            # 无有效数据时，根据查询类型显示对应提示
            query_type = "激素" if hormone_name else "基因" if raw_gene_query else "疾病"
            query_target = hormone_name or raw_gene_query or disease_name
            if query_type == "疾病":
                summary = f"数据库中未匹配到{query_target}相关的激素-基因-受体数据（可能是疾病名称格式不规范或数据关联缺失）"
            else:
                summary = f"数据库中存在{query_target.upper() if query_type == '基因' else query_target}相关记录，但未匹配到关联的激素-受体数据（可能是查询条件过严或数据关联缺失）"
    
    except Exception as e:
        logger.error(f"数据查询失败：{str(e)}", exc_info=True)
        return JsonResponse({
            "code": 500,
            "msg": f"数据获取失败：{str(e)[:100]}"
        })
    
    # 7. 整理模板数据（适配三种查询类型，优化标题显示）
    query_type = "激素" if hormone_name else "基因" if raw_gene_query else "疾病"
    display_target = hormone_name or (raw_gene_query.upper() if raw_gene_query else standardize_disease_name(disease_name))
    template_data = {
        "query_type": query_type,  # 新增：传递查询类型给PDF模板
        "display_target": display_target,  # 新增：传递格式化的查询目标给PDF模板
        "hormone_name": clean_text(hormone_name) or "未知",
        "gene_name": clean_text(raw_gene_query.upper()) if raw_gene_query else "未知",  # 新增：基因名称单独传递
        "disease_name": clean_text(standardize_disease_name(disease_name)) if disease_name else "未知",  # 新增：疾病名称单独传递
        "match_type": "精确匹配" if match_type == "exact" else "模糊匹配",
        "current_time": timezone.now().strftime("%Y-%m-%d %H:%M:%S"),
        "related_genes": real_related_genes if real_related_genes else [
            {"Hormone_Name": "N/A", "Related_Genes": "N/A", "PMID": "N/A", "Related_Diseases": "N/A", "Gene_Sequence": "N/A"}
        ],
        "receptor_full": real_receptor_full if real_receptor_full else [
            {"Hormone_Name": "N/A", "Hormone_Species_Name": "N/A", "Hormone_UniProt_ID": "N/A",
             "Receptor_Name": "N/A", "Receptor_Species_Name": "N/A", "Receptor_Coding_Genes": "N/A"}
        ],
        "receptor_info": real_receptor_info if real_receptor_info else [
            {"Hormone_Name": "N/A", "PubChem_ID": "N/A", "Receptor_Name": "N/A",
             "Receptor_UniProt_ID": "N/A", "Receptor_Species_Name": "N/A", "Receptor_Coding_Genes": "N/A"}
        ],
        "total_records": total_records,
        "summary": summary or "无分析总结数据"
    }
    
    # 8. 生成PDF并下载（优化文件名，适配三种查询类型）
    try:
        pdf = render_to_pdf('hormone_report.html', template_data)
        if pdf:
            # 优化文件名：根据查询类型生成更具辨识度的文件名
            base_name = ""
            if hormone_name:
                base_name = standardize_hormone_name(hormone_name)
            elif raw_gene_query:
                base_name = standardize_gene_name(raw_gene_query)
            elif disease_name:
                base_name = standardize_disease_name(disease_name).replace(" ", "_")
            
            # 安全处理文件名，避免非法字符
            safe_name = re.sub(r'[^\w\s-]', '', base_name).strip().replace(' ', '_') or "Unknown_Query"
            filename = f"{safe_name}_{query_type}_{match_type}_report.pdf"
            
            response = HttpResponse(pdf, content_type='application/pdf')
            response['Content-Disposition'] = (
                f'attachment; filename="{urllib.parse.quote(filename)}"; '
                f'filename*=UTF-8\'\'{urllib.parse.quote(filename)}'
            )
            response['Access-Control-Expose-Headers'] = 'Content-Disposition'
            return response
        else:
            return JsonResponse({"code": 500, "msg": "PDF生成失败：模板渲染为空"})
    except Exception as e:
        logger.error(f"PDF生成失败：{str(e)}", exc_info=True)
        return JsonResponse({"code": 500, "msg": f"PDF生成失败：{str(e)[:100]}"})
        
@csrf_exempt  # 如果已通过 CSRF token 验证，可去掉此装饰器
def extract_entity_for_pdf_view(request):
    """
    接收用户问题，返回提取到的实体列表（用于 PDF 生成）
    """
    if request.method != 'POST':
        return JsonResponse({'code': 405, 'msg': 'Method not allowed'})

    try:
        data = json.loads(request.body)
        question = data.get('question', '').strip()
        if not question:
            return JsonResponse({'code': 400, 'msg': 'Question is required'})

        entities = extract_entities_for_pdf(question)

        return JsonResponse({
            'code': 200,
            'msg': 'success',
            'data': entities
        })
    except json.JSONDecodeError:
        return JsonResponse({'code': 400, 'msg': 'Invalid JSON'})
    except Exception as e:
        logger.error(f"Entity extraction view error: {str(e)}")
        return JsonResponse({'code': 500, 'msg': str(e)})
# ===================== DeepSeek聊天视图 =====================
class UnifiedDeepSeekView(View):
    """
    统一的DeepSeek接口 - 内部自动进行RAG检索增强
    """
    
    def post(self, request):
        request_id = datetime.now().strftime('%Y%m%d%H%M%S%f')
        logger.info(f"[{request_id}] ========== New DeepSeek Request ==========")
        
        try:
            # 解析请求数据
            try:
                request_data = json.loads(request.body) if request.body else {}
                logger.info(f"[{request_id}] Request body parsed: {request_data}")
            except json.JSONDecodeError as e:
                logger.error(f"[{request_id}] JSON decode error: {str(e)}")
                request_data = request.POST.dict()
                logger.info(f"[{request_id}] Using POST data: {request_data}")
            
            question = request_data.get("input_text", "").strip()
            logger.info(f"[{request_id}] Question: '{question}'")
            
            if not question:
                logger.warning(f"[{request_id}] Empty question received")
                return JsonResponse({
                    "code": 400, 
                    "msg": "Please enter a valid question"
                })
            
            # 执行RAG检索增强
            logger.info(f"[{request_id}] Starting RAG search...")
            rag_result = rag_deepseek.rag_search(question)
            
            logger.info(f"[{request_id}] RAG search completed")
            logger.info(f"[{request_id}] Result: success={rag_result['success']}, has_data={rag_result.get('has_data', False)}")
            
            if rag_result["success"]:
                response_data = {
                    "code": 200,
                    "msg": "Success",
                    "data": {
                        "question": question,
                        "answer": rag_result["answer"],
                        "model": rag_result.get("model", "deepseek-v4-pro"),
                        "has_data": rag_result.get("has_data", False),
                        "context_type": rag_result.get("context_type", "unknown"),
                    }
                }
                
                # 如果有token使用信息，也返回
                if rag_result.get("token_usage"):
                    response_data["data"]["token_usage"] = rag_result["token_usage"]
                    logger.info(f"[{request_id}] Token usage: {rag_result['token_usage']}")
                
                logger.info(f"[{request_id}] Response prepared: context_type={response_data['data']['context_type']}")
                logger.info(f"[{request_id}] Answer preview: {response_data['data']['answer'][:200]}...")
                logger.info(f"[{request_id}] ========== Request Completed ==========")
                
                return JsonResponse(response_data)
            else:
                logger.error(f"[{request_id}] RAG search returned success=False")
                return JsonResponse({
                    "code": 500,
                    "msg": "Failed to generate response",
                    "data": {
                        "question": question,
                        "answer": "Sorry, I encountered an error while processing your request.",
                        "has_data": False
                    }
                })
            
        except Exception as e:
            logger.error(f"[{request_id}] System exception: {str(e)}")
            logger.error(traceback.format_exc())
            return JsonResponse({
                "code": 500,
                "msg": f"System error: {str(e)[:100]}"
            })


# ===================== 基因分析工具类 =====================
class GeneAnalysisTool:
    def __init__(self, task):
        self.task = task
        self.expression_data = None
        self.sample_info = None
        self.gene_annotations = None
        self.results_dir = os.path.join(settings.MEDIA_ROOT, 'gene_results', str(task.task_id))
        os.makedirs(self.results_dir, exist_ok=True)
        
        plt.rcParams["font.family"] = ["SimHei", "WenQuanYi Micro Hei", "Heiti TC"]
        plt.rcParams['axes.unicode_minus'] = False
    
    def update_progress(self, progress, status_text=None):
        self.task.progress = progress
        if status_text:
            self.task.status = status_text
        self.task.save()
    
    def load_data(self):
        self.update_progress(10, 'running')
        
        expr_file = self.task.uploaded_files.filter(file_type='expression-data').first()
        if not expr_file:
            raise ValueError("未找到表达数据文件")
        
        file_ext = os.path.splitext(expr_file.filename)[1].lower()
        if file_ext == '.csv':
            self.expression_data = pd.read_csv(expr_file.file.path, index_col=0)
        elif file_ext == '.tsv' or file_ext == '.txt':
            self.expression_data = pd.read_csv(expr_file.file.path, sep='\t', index_col=0)
        elif file_ext == '.xlsx':
            self.expression_data = pd.read_excel(expr_file.file.path, index_col=0)
        else:
            raise ValueError(f"不支持的表达数据格式: {file_ext}")
        
        sample_file = self.task.uploaded_files.filter(file_type='sample-info').first()
        if sample_file:
            file_ext = os.path.splitext(sample_file.filename)[1].lower()
            if file_ext == '.csv':
                self.sample_info = pd.read_csv(sample_file.file.path, index_col=0)
            elif file_ext == '.tsv' or file_ext == '.txt':
                self.sample_info = pd.read_csv(sample_file.file.path, sep='\t', index_col=0)
            elif file_ext == '.xlsx':
                self.sample_info = pd.read_excel(sample_file.file.path, index_col=0)
        
        self.update_progress(20, 'running')
        return True
    
    def normalize_data(self):
        self.update_progress(25, 'running')
        
        method = self.task.normalization_params.get('method', 'tpm')
        
        if method == 'tpm':
            row_sums = self.expression_data.sum(axis=1)
            self.expression_data = self.expression_data.div(row_sums, axis=0) * 1e6
            self.expression_data = np.log2(self.expression_data + 1)
        elif method == 'zscore':
            self.expression_data = (self.expression_data - self.expression_data.mean()) / self.expression_data.std()
        elif method == 'quantile':
            rank_mean = self.expression_data.stack().groupby(self.expression_data.rank(method='first').stack().astype(int)).mean()
            self.expression_data = self.expression_data.rank(method='min').stack().astype(int).map(rank_mean).unstack()
        
        self.update_progress(35, 'running')
        return True
    
    def differential_expression(self):
        self.update_progress(40, 'running')
        
        group_col = self.task.de_params.get('group_column', 'group')
        group1 = self.task.de_params.get('group1', 'control')
        group2 = self.task.de_params.get('group2', 'treatment')
        
        if not self.sample_info or group_col not in self.sample_info.columns:
            raise ValueError(f"样本信息中未找到分组列: {group_col}")
        
        samples1 = self.sample_info[self.sample_info[group_col] == group1].index.tolist()
        samples2 = self.sample_info[self.sample_info[group_col] == group2].index.tolist()
        
        samples1 = [s for s in samples1 if s in self.expression_data.columns]
        samples2 = [s for s in samples2 if s in self.expression_data.columns]
        
        if len(samples1) == 0 or len(samples2) == 0:
            raise ValueError("未找到有效的样本数据用于差异分析")
        
        de_results = []
        for gene, row in self.expression_data.iterrows():
            expr1 = row[samples1].values
            expr2 = row[samples2].values
            
            stat, p_value = stats.ttest_ind(expr1, expr2, equal_var=False)
            
            mean1 = np.mean(expr1)
            mean2 = np.mean(expr2)
            log2fc = mean2 - mean1
            
            de_results.append({
                'gene': gene,
                'mean_group1': mean1,
                'mean_group2': mean2,
                'log2fc': log2fc,
                'p_value': p_value,
                'significant': 'yes' if p_value < 0.05 and abs(log2fc) > 1 else 'no'
            })
        
        self.de_df = pd.DataFrame(de_results)
        self.de_df['adj_p_value'] = stats.false_discovery_control(self.de_df['p_value'])
        
        self.update_progress(60, 'running')
        return True
    
    def generate_visualizations(self):
        self.update_progress(65, 'running')
        
        plt.figure(figsize=(10, 8))
        sig = self.de_df[(self.de_df['adj_p_value'] < 0.05) & (abs(self.de_df['log2fc']) > 1)]
        non_sig = self.de_df[(self.de_df['adj_p_value'] >= 0.05) | (abs(self.de_df['log2fc']) <= 1)]
        
        plt.scatter(non_sig['log2fc'], -np.log10(non_sig['adj_p_value']), 
                   color='gray', alpha=0.5, label='不显著')
        plt.scatter(sig['log2fc'], -np.log10(sig['adj_p_value']), 
                   color='red', alpha=0.7, label='显著')
        
        plt.axhline(y=-np.log10(0.05), color='black', linestyle='--')
        plt.axvline(x=1, color='black', linestyle='--')
        plt.axvline(x=-1, color='black', linestyle='--')
        
        plt.xlabel('log2(倍数变化)')
        plt.ylabel('-log10(校正后p值)')
        plt.title('差异表达基因火山图')
        plt.legend()
        plt.tight_layout()
        
        volcano_path = os.path.join(self.results_dir, 'volcano_plot.png')
        plt.savefig(volcano_path, dpi=300)
        plt.close()
        
        top_genes = self.de_df.sort_values('adj_p_value').head(50)['gene'].tolist()
        heatmap_path = None
        if len(top_genes) >= 5:
            plt.figure(figsize=(12, 10))
            heatmap_data = self.expression_data.loc[top_genes]
            heatmap_data = (heatmap_data - heatmap_data.mean(axis=1)[:, np.newaxis]) / heatmap_data.std(axis=1)[:, np.newaxis]
            
            sns.heatmap(heatmap_data, cmap='coolwarm', center=0)
            plt.title('Top 50 差异表达基因热图')
            plt.tight_layout()
            
            heatmap_path = os.path.join(self.results_dir, 'heatmap.png')
            plt.savefig(heatmap_path, dpi=300)
            plt.close()
        
        self.update_progress(80, 'running')
        return {
            'volcano_plot': volcano_path,
            'heatmap': heatmap_path
        }
    
    def save_results(self, visualizations):
        results = []
        
        norm_path = os.path.join(self.results_dir, 'normalized_expression.csv')
        self.expression_data.to_csv(norm_path)
        
        with open(norm_path, 'rb') as f:
            content = f.read()
            file_name = os.path.basename(norm_path)
            file_path = default_storage.save(f'gene_results/{self.task.task_id}/{file_name}', ContentFile(content))
            results.append(GeneAnalysisResult(
                task=self.task,
                result_type='normalized-data',
                file=file_path,
                filename=file_name
            ))
        
        de_path = os.path.join(self.results_dir, 'differential_expression_results.csv')
        self.de_df.to_csv(de_path, index=False)
        
        with open(de_path, 'rb') as f:
            content = f.read()
            file_name = os.path.basename(de_path)
            file_path = default_storage.save(f'gene_results/{self.task.task_id}/{file_name}', ContentFile(content))
            results.append(GeneAnalysisResult(
                task=self.task,
                result_type='de-results',
                file=file_path,
                filename=file_name
            ))
        
        if 'volcano_plot' in visualizations and visualizations['volcano_plot']:
            with open(visualizations['volcano_plot'], 'rb') as f:
                content = f.read()
                file_name = os.path.basename(visualizations['volcano_plot'])
                file_path = default_storage.save(f'gene_results/{self.task.task_id}/{file_name}', ContentFile(content))
                results.append(GeneAnalysisResult(
                    task=self.task,
                    result_type='volcano-plot',
                    file=file_path,
                    filename=file_name
                ))
        
        if 'heatmap' in visualizations and visualizations['heatmap']:
            with open(visualizations['heatmap'], 'rb') as f:
                content = f.read()
                file_name = os.path.basename(visualizations['heatmap'])
                file_path = default_storage.save(f'gene_results/{self.task.task_id}/{file_name}', ContentFile(content))
                results.append(GeneAnalysisResult(
                    task=self.task,
                    result_type='heatmap',
                    file=file_path,
                    filename=file_name
                ))
        
        GeneAnalysisResult.objects.bulk_create(results)
        
        self.update_progress(100, 'completed')
        self.task.completed_at = datetime.now()
        self.task.save()
        
        return results

# ===================== 错误处理视图 =====================
def bad_request(request, exception=None):
    if request.is_ajax() or request.headers.get('X-Requested-With') == 'XMLHttpRequest':
        return JsonResponse({'error': '错误的请求'}, status=400)
    return render(request, '400.html', status=400)

def permission_denied(request, exception=None):
    if request.is_ajax() or request.headers.get('X-Requested-With') == 'XMLHttpRequest':
        return JsonResponse({'error': '权限不足'}, status=403)
    return render(request, '403.html', status=403)

def page_not_found(request, exception=None):
    if request.is_ajax() or request.headers.get('X-Requested-With') == 'XMLHttpRequest':
        return JsonResponse({'error': '页面未找到'}, status=404)
    return render(request, '404.html', status=404)

def server_error(request):
    if request.is_ajax() or request.headers.get('X-Requested-With') == 'XMLHttpRequest':
        return JsonResponse({'error': '服务器内部错误'}, status=500)
    return render(request, '500.html', status=500)

# ===================== 注册、登录、登出相关视图 =====================
@csrf_protect
def login_view(request):
    if request.method == 'POST':
        username = request.POST.get('username', '').strip()
        password = request.POST.get('password', '').strip()

        if not username or not password:
            return JsonResponse(
                {'success': False, 'error': '用户名和密码不能为空'},
                status=400
            )

        fail_key = f"login_fail_{username}"
        fail_count = cache.get(fail_key, 0)
        if fail_count >= 5:
            return JsonResponse(
                {'success': False, 'error': '多次登录失败，账号已临时锁定5分钟'},
                status=429
            )

        user = authenticate(request, username=username, password=password)

        if user:
            if user.is_active:
                login(request, user)
                cache.delete(fail_key)
                logger.info(f"用户 {username} 登录成功")
                return JsonResponse({
                    'success': True,
                    'username': user.username,
                    'is_admin': user.is_staff,
                })
            else:
                logger.warning(f"禁用用户 {username} 尝试登录")
                return JsonResponse(
                    {'success': False, 'error': '账号已被禁用'},
                    status=403
                )
        else:
            cache.set(fail_key, fail_count + 1, 300)
            logger.warning(f"用户 {username} 登录失败（剩余尝试次数：{5 - fail_count - 1}）")
            return JsonResponse(
                {'success': False, 'error': '用户名或密码错误'},
                status=401
            )

    return render(request, 'login_view.html')

@csrf_protect
def register_view(request):
    if request.method == 'POST':
        username = request.POST.get('username', '').strip()
        password = request.POST.get('password', '').strip()

        errors = []
        if not username or not password:
            errors.append('用户名和密码不能为空')
        if User.objects.filter(username=username).exists():
            errors.append('用户名已存在')
        if len(password) < 6:
            errors.append('密码长度至少6位')
        if not re.search(r'[A-Za-z]', password) or not re.search(r'[0-9]', password):
            errors.append('密码需包含字母和数字')
        if not username.isalnum():
            errors.append('用户名仅支持字母和数字')

        if errors:
            return JsonResponse(
                {'success': False, 'error': '; '.join(errors)},
                status=400
            )

        try:
            user = User.objects.create_user(username=username, password=password)
            logger.info(f"用户 {username} 注册成功")
            return JsonResponse({'success': True})
        except Exception as e:
            logger.error(f"注册失败: {str(e)}")
            return JsonResponse({'success': False, 'error': '注册失败，请稍后重试'})
    else:
        return render(request, 'register.html')

def logout_view(request):
    logout(request)
    return redirect('hormone_app:login_view')

# ===================== 页面视图 =====================
def page1(request):
    return render(request, '1.html')

def page2(request):
    return render(request, '2.html')

def page3(request):
    return render(request, '3.html')

def page4(request):
    return render(request, '4.html')

def page5(request):
    return render(request, '5.html')

def page6(request):
    return render(request, '6.html')

# ===================== 搜索功能视图 =====================
def _safe_text(value, fallback='N/A'):
    if value is None:
        return fallback
    value = str(value).strip()
    return value if value else fallback


def _truncate_text(value, length=240):
    value = _safe_text(value, '')
    return value[:length] + ('...' if len(value) > length else '')


def _entity_url(entity_type, key):
    return f"/hormone_app/entity/{entity_type}/?key={quote(str(key or ''))}"


def _append_entity(results, seen, entity_type, title, subtitle='', description='', key=None, count=0, source=''):
    title = _safe_text(title, '')
    key = _safe_text(key if key is not None else title, '')
    if not title or not key:
        return
    dedupe_key = (entity_type, key.lower())
    if dedupe_key in seen:
        return
    seen.add(dedupe_key)
    badge_map = {'hormone': 'Hormone', 'gene': 'Gene', 'disease': 'Disease', 'receptor': 'Receptor'}
    results.append({
        'type': entity_type,
        'type_label': badge_map.get(entity_type, entity_type.title()),
        'title': title,
        'subtitle': subtitle,
        'description': _truncate_text(description),
        'count': count,
        'source': source,
        'url': _entity_url(entity_type, key),
    })


def search(request):
    """NCBI-like entity search across hormones, genes, diseases and receptors."""
    query = request.GET.get('query', '').strip()
    entity_filter = request.GET.get('type', 'all').strip().lower() or 'all'
    allowed_types = {'all', 'hormone', 'gene', 'disease', 'receptor'}
    if entity_filter not in allowed_types:
        entity_filter = 'all'

    results = []
    seen = set()
    counts = {'hormone': 0, 'gene': 0, 'disease': 0, 'receptor': 0}

    if query:
        if entity_filter in ('all', 'hormone'):
            hormone_rows = HormoneRelatedGene.objects.filter(hormone_name__icontains=query).exclude(hormone_name__isnull=True).exclude(hormone_name='').values('hormone_name').annotate(record_count=Count('id')).order_by('-record_count', 'hormone_name')[:80]
            for row in hormone_rows:
                name = row['hormone_name']
                sample = HormoneRelatedGene.objects.filter(hormone_name=name).first()
                _append_entity(results, seen, 'hormone', name, subtitle=f"{row['record_count']} hormone-gene-disease records", description=f"Organism: {_safe_text(getattr(sample, 'organism', ''))}; example gene: {_safe_text(getattr(sample, 'related_genes', ''))}; example disease: {_safe_text(getattr(sample, 'related_diseases', ''))}", count=row['record_count'], source='hormone_data')
            for model, source in ((HormoneReceptorFull, 'peptide receptor table'), (HormoneReceptorInfo, 'non-peptide receptor table')):
                for row in model.objects.filter(hormone_name__icontains=query).exclude(hormone_name__isnull=True).exclude(hormone_name='').values('hormone_name').annotate(record_count=Count('id')).order_by('-record_count', 'hormone_name')[:30]:
                    _append_entity(results, seen, 'hormone', row['hormone_name'], subtitle=f"{row['record_count']} receptor records", description=source, count=row['record_count'], source=source)

        if entity_filter in ('all', 'gene'):
            gene_rows = HormoneRelatedGene.objects.filter(related_genes__icontains=query).exclude(related_genes__isnull=True).exclude(related_genes='').values('related_genes').annotate(record_count=Count('id')).order_by('-record_count', 'related_genes')[:80]
            for row in gene_rows:
                gene = row['related_genes']
                sample = HormoneRelatedGene.objects.filter(related_genes=gene).first()
                _append_entity(results, seen, 'gene', gene, subtitle=f"{row['record_count']} related hormone/disease records", description=f"Example hormone: {_safe_text(getattr(sample, 'hormone_name', ''))}; chromosome: {_safe_text(getattr(sample, 'chromosome', ''))}; description: {_safe_text(getattr(sample, 'gene_description', ''))}", count=row['record_count'], source='hormone_data')

        if entity_filter in ('all', 'disease'):
            disease_rows = HormoneRelatedGene.objects.filter(Q(related_diseases__icontains=query) | Q(doid_standardized_name__icontains=query) | Q(DO_ID__icontains=query)).exclude(related_diseases__isnull=True).exclude(related_diseases='').values('related_diseases').annotate(record_count=Count('id')).order_by('-record_count', 'related_diseases')[:80]
            for row in disease_rows:
                disease = row['related_diseases']
                sample = HormoneRelatedGene.objects.filter(related_diseases=disease).first()
                _append_entity(results, seen, 'disease', disease, subtitle=f"{row['record_count']} gene-disease records", description=f"DOID: {_safe_text(getattr(sample, 'DO_ID', ''))}; standardized name: {_safe_text(getattr(sample, 'doid_standardized_name', ''))}; example gene: {_safe_text(getattr(sample, 'related_genes', ''))}", count=row['record_count'], source='hormone_data')

        if entity_filter in ('all', 'receptor'):
            receptor_sources = [(Hormone, 'receptor_name', 'classic receptor table'), (HormoneReceptorFull, 'receptor_name', 'peptide receptor table'), (HormoneReceptorInfo, 'receptor_name', 'non-peptide receptor table')]
            for model, field, source in receptor_sources:
                q = Q(**{f'{field}__icontains': query})
                if model is Hormone:
                    q |= Q(gene_name_symbol__icontains=query) | Q(ligand__icontains=query)
                else:
                    q |= Q(receptor_uniprot_id__icontains=query) | Q(receptor_coding_genes__icontains=query)
                for row in model.objects.filter(q).exclude(**{f'{field}__isnull': True}).exclude(**{field: ''}).values(field).annotate(record_count=Count('id')).order_by('-record_count', field)[:50]:
                    name = row[field]
                    sample = model.objects.filter(**{field: name}).first()
                    if model is Hormone:
                        desc = f"Gene symbol: {_safe_text(getattr(sample, 'gene_name_symbol', ''))}; ligand: {_safe_text(getattr(sample, 'ligand', ''))}; family: {_safe_text(getattr(sample, 'gene_family', ''))}"
                    else:
                        desc = f"UniProt: {_safe_text(getattr(sample, 'receptor_uniprot_id', ''))}; gene: {_safe_text(getattr(sample, 'receptor_coding_genes', ''))}; species: {_safe_text(getattr(sample, 'receptor_species_name', ''))}"
                    _append_entity(results, seen, 'receptor', name, subtitle=f"{row['record_count']} receptor records", description=desc, count=row['record_count'], source=source)

    for item in results:
        counts[item['type']] += 1

    return render(request, 'search_results.html', {'query': query, 'entity_filter': entity_filter, 'results': results, 'total': len(results), 'counts': counts})


def entity_detail(request, entity_type):
    """Entity detail page for a hormone/gene/disease/receptor search hit."""
    entity_type = (entity_type or '').lower()
    key = request.GET.get('key', '').strip()
    if entity_type not in {'hormone', 'gene', 'disease', 'receptor'} or not key:
        return render(request, 'entity_detail.html', {'error': 'Invalid entity request', 'entity_type': entity_type, 'key': key})

    context = {'entity_type': entity_type, 'entity_label': entity_type.title(), 'key': key, 'records': [], 'receptor_records': [], 'stats': {}}

    if entity_type == 'hormone':
        records = HormoneRelatedGene.objects.filter(hormone_name=key).order_by('id')[:200]
        receptor_records = list(HormoneReceptorFull.objects.filter(hormone_name=key)[:100]) + list(HormoneReceptorInfo.objects.filter(hormone_name=key)[:100])
        all_records = HormoneRelatedGene.objects.filter(hormone_name=key)
        sample = all_records.first()
        entity_info = {}
        if sample:
            entity_info = {
                'hormone_name': getattr(sample, 'hormone_name', ''),
                'hormone_accession': getattr(sample, 'hormone_accession', ''),
                'organism': getattr(sample, 'organism', ''),
            }
        pep = HormoneReceptorFull.objects.filter(hormone_name=key).first()
        if pep:
            entity_info.update({
                'hormone_uniprot_id': getattr(pep, 'hormone_uniprot_id', ''),
                'hormone_species_name': getattr(pep, 'hormone_species_name', ''),
                'hormone_coding_genes': getattr(pep, 'hormone_coding_genes', ''),
                'hormone_coding_genes_sequence': getattr(pep, 'hormone_coding_genes_sequence', ''),
            })
        linked_genes = list(all_records.exclude(related_genes__isnull=True).exclude(related_genes='').values_list('related_genes', flat=True).distinct().order_by('related_genes')[:500])
        linked_diseases = list(all_records.exclude(related_diseases__isnull=True).exclude(related_diseases='').values_list('related_diseases', flat=True).distinct().order_by('related_diseases')[:500])
        linked_receptors = list(set(
            [r.receptor_name for r in HormoneReceptorFull.objects.filter(hormone_name=key).exclude(receptor_name__isnull=True).exclude(receptor_name='')[:100]] +
            [r.receptor_name for r in HormoneReceptorInfo.objects.filter(hormone_name=key).exclude(receptor_name__isnull=True).exclude(receptor_name='')[:100]]
        ))
        context.update({'records': records, 'receptor_records': receptor_records, 'entity_info': entity_info, 'linked_genes': linked_genes, 'linked_diseases': linked_diseases, 'linked_receptors': linked_receptors, 'stats': {'Records': all_records.count(), 'Genes': len(linked_genes), 'Diseases': len(linked_diseases), 'Organisms': all_records.values('organism').distinct().count()}})
    elif entity_type == 'gene':
        records = HormoneRelatedGene.objects.filter(related_genes=key).order_by('id')[:200]
        all_records = HormoneRelatedGene.objects.filter(related_genes=key)
        # Build NCBI-style gene info card from first available record
        sample = all_records.first()
        gene_info = {}
        if sample:
            gene_info = {
                'gene_name': getattr(sample, 'related_genes', ''),
                'description': getattr(sample, 'gene_description', ''),
                'chromosome': getattr(sample, 'chromosome', ''),
                'map_location': getattr(sample, 'map_location', ''),
                'synonyms': getattr(sample, 'gene_synonyms', ''),
                'chr_accession': getattr(sample, 'chr_accession', ''),
                'hormone_accession': getattr(sample, 'hormone_accession', ''),
                'protein_accession': getattr(sample, 'protein_accession', ''),
                'protein_sequence': getattr(sample, 'protein_sequence', ''),
                'mrna_accession': getattr(sample, 'mrna_accession', ''),
                'mrna_sequence': getattr(sample, 'mrna_sequence', ''),
                'organism': getattr(sample, 'organism', ''),
                'evidence_source': getattr(sample, 'evidence_source', ''),
                'regulation_type': getattr(sample, 'regulation_type', ''),
                'relationship_gene_disease': getattr(sample, 'relationship_gene_disease', ''),
            }
        linked_hormones = list(all_records.exclude(hormone_name__isnull=True).exclude(hormone_name='').values_list('hormone_name', flat=True).distinct().order_by('hormone_name')[:500])
        linked_diseases = list(all_records.exclude(related_diseases__isnull=True).exclude(related_diseases='').values_list('related_diseases', flat=True).distinct().order_by('related_diseases')[:500])
        context.update({
            'records': records,
            'gene_info': gene_info,
            'linked_hormones': linked_hormones,
            'linked_diseases': linked_diseases,
            'stats': {
                'Records': all_records.count(),
                'Hormones': len(linked_hormones),
                'Diseases': len(linked_diseases),
                'Organisms': all_records.values('organism').distinct().count(),
            }
        })
    elif entity_type == 'disease':
        records = HormoneRelatedGene.objects.filter(related_diseases=key).order_by('id')[:200]
        all_records = HormoneRelatedGene.objects.filter(related_diseases=key)
        sample = all_records.first()
        entity_info = {}
        if sample:
            entity_info = {
                'disease_name': getattr(sample, 'related_diseases', ''),
                'DO_ID': getattr(sample, 'DO_ID', ''),
                'doid_standardized_name': getattr(sample, 'doid_standardized_name', ''),
                'DO_Standardized_Terminology': getattr(sample, 'DO_Standardized_Terminology', ''),
                'DO_Match_disease_names': getattr(sample, 'DO_Match_disease_names', ''),
                'relationship_gene_disease': getattr(sample, 'relationship_gene_disease', ''),
                'has_pmid_gene_disease': getattr(sample, 'has_pmid_gene_disease', ''),
                'evidence_source': getattr(sample, 'evidence_source', ''),
            }
        linked_hormones = list(all_records.exclude(hormone_name__isnull=True).exclude(hormone_name='').values_list('hormone_name', flat=True).distinct().order_by('hormone_name')[:500])
        linked_genes = list(all_records.exclude(related_genes__isnull=True).exclude(related_genes='').values_list('related_genes', flat=True).distinct().order_by('related_genes')[:500])
        context.update({'records': records, 'entity_info': entity_info, 'linked_hormones': linked_hormones, 'linked_genes': linked_genes, 'stats': {'Records': all_records.count(), 'Hormones': len(linked_hormones), 'Genes': len(linked_genes), 'DO Terms': all_records.values('DO_ID').distinct().count()}})
    else:
        classic = Hormone.objects.filter(receptor_name=key)[:100]
        peptide = HormoneReceptorFull.objects.filter(receptor_name=key)[:100]
        non_peptide = HormoneReceptorInfo.objects.filter(receptor_name=key)[:100]
        receptor_records = list(classic) + list(peptide) + list(non_peptide)
        entity_info = {}
        cls = Hormone.objects.filter(receptor_name=key).first()
        if cls:
            entity_info.update({'receptor_type': getattr(cls, 'receptor_type', ''), 'receptor_name': getattr(cls, 'receptor_name', ''), 'gene_family': getattr(cls, 'gene_family', ''), 'gene_name_symbol': getattr(cls, 'gene_name_symbol', ''), 'ligand': getattr(cls, 'ligand', ''), 'ligand_type_comments': getattr(cls, 'ligand_type_comments', '')})
        np = HormoneReceptorInfo.objects.filter(receptor_name=key).first()
        if np:
            entity_info.update({'pubchem_id': getattr(np, 'pubchem_id', ''), 'receptor_uniprot_id': getattr(np, 'receptor_uniprot_id', ''), 'receptor_coding_genes': getattr(np, 'receptor_coding_genes', ''), 'receptor_species_name': getattr(np, 'receptor_species_name', ''), 'receptor_coding_genes_sequence': getattr(np, 'receptor_coding_genes_sequence', '')})
        linked_hormones = list(set(
            [r.hormone_name for r in HormoneReceptorFull.objects.filter(receptor_name=key).exclude(hormone_name__isnull=True).exclude(hormone_name='')[:200]] +
            [r.hormone_name for r in HormoneReceptorInfo.objects.filter(receptor_name=key).exclude(hormone_name__isnull=True).exclude(hormone_name='')[:200]]
        ))
        linked_genes = list(set(
            [r.receptor_coding_genes for r in HormoneReceptorInfo.objects.filter(receptor_name=key).exclude(receptor_coding_genes__isnull=True).exclude(receptor_coding_genes='')[:200]] +
            [r.receptor_coding_genes for r in HormoneReceptorFull.objects.filter(receptor_name=key).exclude(receptor_coding_genes__isnull=True).exclude(receptor_coding_genes='')[:200]]
        ))
        context.update({'receptor_records': receptor_records, 'entity_info': entity_info, 'linked_hormones': sorted(linked_hormones), 'linked_genes': sorted(linked_genes), 'stats': {'Records': len(receptor_records), 'Classic': Hormone.objects.filter(receptor_name=key).count(), 'Peptide': HormoneReceptorFull.objects.filter(receptor_name=key).count(), 'Non-peptide': HormoneReceptorInfo.objects.filter(receptor_name=key).count()}})

    return render(request, 'entity_detail.html', context)

# ===================== 基础视图 =====================
def hormone_receptor_list(request):
    try:
        receptors = Hormone.objects.all()
        return render(
            request,
            'receptors/list.html',
            {'receptors': receptors}
        )
    except Exception as e:
        logger.error(f"获取激素受体列表失败: {str(e)}", exc_info=True)
        if request.headers.get('X-Requested-With') == 'XMLHttpRequest':
            return JsonResponse({'success': False, 'error': '数据加载失败'}, status=500)
        return render(request, 'error.html', {'error': '系统繁忙，请稍后重试'})

def home(request):
    """
    首页视图 - 添加向量存储初始化
    """
    try:
        # 注意：这里直接使用已导入的vector_store和initialize_vector_store_async
        # 不再重新导入，避免循环导入问题
        
        # 检查是否需要初始化向量存储
        # 条件：未加载、未在加载中、并且模型可用
        should_initialize = (
            not vector_store.is_loaded and 
            not vector_store.loading and
            hasattr(vector_store, 'model_name')
        )
        
        if should_initialize:
            try:
                # 启动后台线程加载向量索引
                initialize_vector_store_async()
                # 修改：使用vector_store的实际模型名称打日志
                actual_model_name = getattr(vector_store, 'model_name', 'unknown-model')
                logger.info(f"语义检索向量索引后台加载已启动（模型：{actual_model_name}）")
            except Exception as e:
                logger.error(f"启动向量索引加载失败: {e}")
        
        # 准备向量存储状态信息（用于调试或前端显示）
        vector_status = {
            "is_loaded": vector_store.is_loaded,
            "is_loading": vector_store.loading,
            # 修改：从vector_store获取实际模型名称，不再硬编码
            "model_name": getattr(vector_store, 'model_name', 'unknown-model'),
            "text_count": len(vector_store.texts) if hasattr(vector_store, 'texts') else 0
        }
        
        # 如果向量索引已加载，记录信息
        if vector_store.is_loaded:
            logger.debug(f"语义检索已就绪，包含 {vector_status['text_count']} 条向量记录")
        
        return render(request, 'home.html', {
            'vector_status': vector_status  # 可选：传递状态到模板
        })
        
    except Exception as e:
        logger.error(f"首页加载失败: {str(e)}", exc_info=True)
        # 出错时仍返回首页，但不包含向量状态
        return render(request, 'home.html')

def list_view(request):
    return render(request, 'list.html')

def docs(request):
    try:
        literature = Literature.objects.filter(
            file__isnull=False
        ).order_by('-upload_date')
        paginator = Paginator(literature, 15)
        page = request.GET.get('page', 1)
        literature = paginator.get_page(page)
        return render(
            request,
            'docs.html',
            {'literature': literature, 'paginator': paginator}
        )
    except Exception as e:
        logger.error(f"获取文献列表失败: {str(e)}")
        return render(request, 'error.html', {'error': '文献数据加载失败'})

def help_view(request):
    return render(request, 'help.html')

def ai_view(request):
    return render(request, 'ai.html')

# ===================== 数据可视化视图 =====================
def data(request):
    try:
        receptors = Hormone.objects.all()
        nodes = []
        links = []
        categories = []

        unique_families = set()
        for receptor in receptors:
            if receptor.gene_family:
                unique_families.add(receptor.gene_family)
        categories = [{'name': family} for family in unique_families]
        categories.append({'name': '配体'})

        node_map = {}
        for receptor in receptors:
            receptor_id = str(receptor.id)
            if receptor_id not in node_map:
                node_map[receptor_id] = {
                    'id': receptor_id,
                    'name': receptor.receptor_name or '未知受体',
                    'symbolSize': 50,
                    'category': receptor.gene_family or '未分类',
                    'value': receptor.receptor_type or ''
                }
                nodes.append(node_map[receptor_id])

            if receptor.ligand:
                ligand_id = f"ligand_{receptor.id}"
                if ligand_id not in node_map:
                    node_map[ligand_id] = {
                        'id': ligand_id,
                        'name': receptor.ligand,
                        'symbolSize': 30,
                        'category': '配体',
                        'value': receptor.ligand_type_comments or '未知'
                    }
                    nodes.append(node_map[ligand_id])
                links.append({
                    'source': receptor_id,
                    'target': ligand_id,
                    'value': '结合'
                })

        graph_data = json.dumps(
            {'nodes': nodes, 'links': links, 'categories': categories},
            ensure_ascii=False
        )

        return render(
            request,
            'data.html',
            {'graph_data': graph_data}
        )

    except Exception as e:
        logger.error(f"生成网络图数据失败: {str(e)}")
        return render(request, 'error.html', {'error': '可视化数据加载失败'})

# ===================== 管理员功能 =====================
def is_admin(user):
    return user.is_authenticated and user.is_staff

@login_required
@user_passes_test(is_admin)
def admin_dashboard(request):
    page = request.GET.get('page', 1)
    page_size = 10

    receptor_full_queryset = HormoneReceptorFull.objects.all().order_by('-id')
    receptor_full_paginator = Paginator(receptor_full_queryset, page_size)
    try:
        receptor_full_list = receptor_full_paginator.page(page)
    except PageNotAnInteger:
        receptor_full_list = receptor_full_paginator.page(1)
    except EmptyPage:
        receptor_full_list = receptor_full_paginator.page(receptor_full_paginator.num_pages)

    receptor_info_queryset = HormoneReceptorInfo.objects.all().order_by('-id')
    receptor_info_paginator = Paginator(receptor_info_queryset, page_size)
    try:
        receptor_info_list = receptor_info_paginator.page(page)
    except PageNotAnInteger:
        receptor_info_list = receptor_info_paginator.page(1)
    except EmptyPage:
        receptor_info_list = receptor_info_paginator.page(receptor_info_paginator.num_pages)

    related_gene_queryset = HormoneRelatedGene.objects.all().order_by('-id')
    related_gene_paginator = Paginator(related_gene_queryset, page_size)
    try:
        related_gene_list = related_gene_paginator.page(page)
    except PageNotAnInteger:
        related_gene_list = related_gene_paginator.page(1)
    except EmptyPage:
        related_gene_list = related_gene_paginator.page(related_gene_paginator.num_pages)

    user_queryset = User.objects.all().order_by('-date_joined')
    user_paginator = Paginator(user_queryset, page_size)
    try:
        users = user_paginator.page(page)
    except PageNotAnInteger:
        users = user_paginator.page(1)
    except EmptyPage:
        users = user_paginator.page(user_paginator.num_pages)

    receptor_full_count = receptor_full_queryset.count()
    receptor_info_count = receptor_info_queryset.count()
    related_gene_count = related_gene_queryset.count()

    user_count = user_queryset.count()
    active_user_count = user_queryset.filter(is_active=True).count()

    today = timezone.now().date()
    today_visits = 156

    context = {
        'receptor_full_list': receptor_full_list,
        'receptor_info_list': receptor_info_list,
        'related_gene_list': related_gene_list,
        'receptor_full_count': receptor_full_count,
        'receptor_info_count': receptor_info_count,
        'related_gene_count': related_gene_count,
        'user_count': user_count,
        'active_user_count': active_user_count,
        'today_visits': today_visits,
        'users': users,
        'current_page': page,
    }

    return render(request, 'admin_dashboard.html', context)

@login_required(login_url='hormone_app:login_view')
def add_user(request):
    if not request.user.is_superuser:
        return JsonResponse(
            {'success': False, 'error': '无权限操作'},
            status=403
        )

    if request.method == 'POST':
        username = request.POST.get('username', '').strip()
        email = request.POST.get('email', '').strip()
        password = request.POST.get('password', '').strip()

        errors = []
        if not username.isalnum():
            errors.append('用户名仅支持字母和数字')
        if not email or '@' not in email:
            errors.append('邮箱格式不正确')
        if User.objects.filter(username=username).exists():
            errors.append('用户名已存在')
        if User.objects.filter(email=email).exists():
            errors.append('邮箱已被使用')
        if len(password) < 6:
            errors.append('密码长度至少6位')

        if errors:
            return JsonResponse(
                {'success': False, 'error': '; '.join(errors)},
                status=400
            )

        try:
            User.objects.create_user(
                username=username,
                email=email,
                password=password
            )
            logger.info(f"管理员 {request.user.username} 创建用户 {username}")
            return JsonResponse({'success': True})
        except Exception as e:
            logger.error(f"创建用户失败: {str(e)}")
            return JsonResponse({'success': False, 'error': '创建失败，请重试'})

    return JsonResponse(
        {'success': False, 'error': '仅支持POST请求'},
        status=405
    )

@login_required(login_url='hormone_app:login_view')
def delete_user(request, user_id):
    if not request.user.is_superuser:
        return JsonResponse(
            {'success': False, 'error': '无权限操作'},
            status=403
        )

    if request.method != 'POST':
        return JsonResponse(
            {'success': False, 'error': '仅支持POST删除'},
            status=405
        )

    try:
        user = User.objects.get(id=user_id)
        if user.is_superuser:
            return JsonResponse(
                {'success': False, 'error': '禁止删除超级管理员'},
                status=403
            )

        username = user.username
        user.delete()
        logger.info(f"管理员 {request.user.username} 删除用户 {username}")
        return JsonResponse({'success': True})
    except User.DoesNotExist:
        return JsonResponse(
            {'success': False, 'error': '用户不存在'},
            status=404
        )
    except Exception as e:
        logger.error(f"删除用户失败: {str(e)}")
        return JsonResponse({'success': False, 'error': '删除失败，请重试'})

@login_required(login_url='hormone_app:login_view')
def upload_literature(request):
    if not request.user.is_staff:
        return JsonResponse(
            {'success': False, 'error': '无上传权限'},
            status=403
        )

    if request.method != 'POST':
        return JsonResponse(
            {'success': False, 'error': '仅支持POST上传'},
            status=405
        )

    title = request.POST.get('title', '').strip()
    authors = request.POST.get('authors', '').strip()
    year = request.POST.get('year', '')
    lit_type = request.POST.get('lit_type', '').strip()
    abstract = request.POST.get('abstract', '').strip()
    file = request.FILES.get('file')

    errors = []
    if not title:
        errors.append('标题不能为空')
    if not authors:
        errors.append('作者不能为空')
    if not year.isdigit() or len(year) != 4:
        errors.append('年份格式错误（需4位数字）')
    if not lit_type:
        errors.append('文献类型不能为空')
    if not file:
        errors.append('请上传文件')
    if file and not file.name.endswith(('.pdf', '.docx')):
        errors.append('仅支持 PDF 或 DOCX 文件')
    if file and file.size > 10 * 1024 * 1024:
        errors.append('文件大小不能超过10MB')

    if errors:
        return JsonResponse(
            {'success': False, 'error': '; '.join(errors)},
            status=400
        )

    try:
        literature = Literature.objects.create(
            title=title,
            authors=authors,
            year=year,
            lit_type=lit_type,
            abstract=abstract,
            file=file,
            upload_date=timezone.now()
        )
        logger.info(f"管理员 {request.user.username} 上传文献: {title}")
        return JsonResponse(
            {'success': True, 'message': '文献上传成功'},
            status=201
        )
    except Exception as e:
        logger.error(f"上传文献失败: {str(e)}")
        return JsonResponse(
            {'success': False, 'error': '上传失败，请重试'},
            status=500
        )

# ===================== 基因分析相关视图 =====================
@login_required
def gene_analysis(request):
    if request.method == 'POST':
        sequence_file = request.FILES.get('sequence_file')
        sequence_text = request.POST.get('sequence_text', '').strip()
        
        if not sequence_file and not sequence_text:
            return JsonResponse(
                {'success': False, 'error': '请上传序列文件或输入序列文本'},
                status=400
            )
        
        try:
            if sequence_file:
                sequence = sequence_file.read().decode('utf-8').strip()
            else:
                sequence = sequence_text
            
            if not validate_sequence(sequence):
                return JsonResponse(
                    {'success': False, 'error': '无效的DNA序列，请检查输入'},
                    status=400
                )
            
            analysis_results = perform_analysis(sequence)
            save_analysis(request.user, sequence, analysis_results)
            
            return JsonResponse({
                'success': True,
                'results': analysis_results
            })
            
        except Exception as e:
            logger.error(f"基因分析失败: {str(e)}", exc_info=True)
            return JsonResponse(
                {'success': False, 'error': '分析过程出错，请重试'},
                status=500
            )
    
    return render(request, 'gene_analysis.html')

@login_required
def analysis_history(request):
    try:
        history = AnalysisHistory.objects.filter(
            user=request.user
        ).order_by('-created_at')
        
        paginator = Paginator(history, 10)
        page = request.GET.get('page', 1)
        history = paginator.get_page(page)
        
        return render(
            request,
            'analysis_history.html',
            {'history': history, 'paginator': paginator}
        )
        
    except Exception as e:
        logger.error(f"获取分析历史失败: {str(e)}")
        return render(
            request,
            'error.html',
            {'error': '无法加载分析历史记录'}
        )

@login_required
def analysis_detail(request, analysis_id):
    try:
        analysis = AnalysisHistory.objects.get(
            id=analysis_id,
            user=request.user
        )
        
        results = json.loads(analysis.results)
        
        return render(
            request,
            'analysis_detail.html',
            {'analysis': analysis, 'results': results}
        )
        
    except AnalysisHistory.DoesNotExist:
        return render(
            request,
            'error.html',
            {'error': '分析记录不存在或您无权访问'}
        )
    except Exception as e:
        logger.error(f"获取分析详情失败: {str(e)}")
        return render(
            request,
            'error.html',
            {'error': '无法加载分析详情'}
        )

def gene_analysis_page(request):
    return render(request, 'hormone_app/gene_analysis.html')

def create_gene_task(request):
    if request.method == 'POST':
        task_name = request.POST.get('task_name', '未命名基因分析')
        task = GeneAnalysisTask.objects.create(
            user=request.user if request.user.is_authenticated else None,
            task_name=task_name
        )
        return JsonResponse({
            'task_id': str(task.task_id),
            'status': task.status,
            'created_at': task.created_at.isoformat()
        })
    return JsonResponse({'error': '无效请求'}, status=400)

def upload_gene_file(request, task_id):
    try:
        task = GeneAnalysisTask.objects.get(task_id=task_id)
    except GeneAnalysisTask.DoesNotExist:
        return JsonResponse({'error': '任务不存在'}, status=404)
    
    if request.method == 'POST' and request.FILES.get('file'):
        file = request.FILES['file']
        file_type = request.POST.get('file_type', 'expression-data')
        
        file_path = default_storage.save(f'gene_uploads/{task_id}/{file.name}', ContentFile(file.read()))
        
        uploaded_file = GeneUploadedFile.objects.create(
            task=task,
            file_type=file_type,
            file=file_path,
            filename=file.name,
            file_size=file.size
        )
        
        return JsonResponse({
            'file_id': str(uploaded_file.file_id),
            'filename': uploaded_file.filename,
            'uploaded_at': uploaded_file.uploaded_at.isoformat()
        })
    
    return JsonResponse({'error': '无效请求'}, status=400)

def run_gene_analysis(request, task_id):
    try:
        task = GeneAnalysisTask.objects.get(task_id=task_id)
    except GeneAnalysisTask.DoesNotExist:
        return JsonResponse({'error': '任务不存在'}, status=404)
    
    if request.method == 'POST':
        task.status = 'running'
        task.progress = 0
        task.save()
        
        try:
            params = json.loads(request.body)
            
            if 'data_params' in params:
                task.data_params = params['data_params']
            if 'normalization_params' in params:
                task.normalization_params = params['normalization_params']
            if 'de_params' in params:
                task.de_params = params['de_params']
            if 'visualization_params' in params:
                task.visualization_params = params['visualization_params']
            task.save()
            
            analysis_tool = GeneAnalysisTool(task)
            analysis_tool.load_data()
            analysis_tool.normalize_data()
            analysis_tool.differential_expression()
            visualizations = analysis_tool.generate_visualizations()
            analysis_tool.save_results(visualizations)
            
            return JsonResponse({
                'status': 'success',
                'task_id': str(task.task_id),
                'progress': task.progress
            })
            
        except Exception as e:
            task.status = 'failed'
            task.error_message = str(e)
            task.save()
            return JsonResponse({
                'error': str(e),
                'status': 'failed'
            }, status=500)
    
    return JsonResponse({'error': '无效请求'}, status=400)

def get_gene_task_status(request, task_id):
    try:
        task = GeneAnalysisTask.objects.get(task_id=task_id)
        return JsonResponse({
            'task_id': str(task.task_id),
            'status': task.status,
            'progress': task.progress,
            'error_message': task.error_message,
            'completed_at': task.completed_at.isoformat() if task.completed_at else None
        })
    except GeneAnalysisTask.DoesNotExist:
        return JsonResponse({'error': '任务不存在'}, status=404)

def get_gene_task_results(request, task_id):
    try:
        task = GeneAnalysisTask.objects.get(task_id=task_id)
        results = task.results.all()
        
        return JsonResponse({
            'task_id': str(task.task_id),
            'status': task.status,
            'results': [{
                'result_id': str(res.result_id),
                'result_type': res.result_type,
                'filename': res.filename,
                'download_url': f'/gene-analysis/download/{res.result_id}/',
                'created_at': res.created_at.isoformat()
            } for res in results]
        })
    except GeneAnalysisTask.DoesNotExist:
        return JsonResponse({'error': '任务不存在'}, status=404)

def download_gene_result(request, result_id):
    try:
        result = GeneAnalysisResult.objects.get(result_id=result_id)
        
        if not os.path.exists(result.file.path):
            return JsonResponse({'error': '文件不存在'}, status=404)
        
        with open(result.file.path, 'rb') as f:
            response = HttpResponse(f.read(), content_type='application/octet-stream')
            response['Content-Disposition'] = f'attachment; filename="{result.filename}"'
            return response
            
    except GeneAnalysisResult.DoesNotExist:
        return JsonResponse({'error': '结果不存在'}, status=404)

def download_all_gene_results(request, task_id):
    try:
        task = GeneAnalysisTask.objects.get(task_id=task_id)
        results = task.results.all()
        
        if not results:
            return JsonResponse({'error': '没有结果文件'}, status=404)
        
        with tempfile.NamedTemporaryFile(delete=False, suffix='.zip') as tmp_file:
            with zipfile.ZipFile(tmp_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
                for res in results:
                    if os.path.exists(res.file.path):
                        zipf.write(res.file.path, res.filename)
            
            with open(tmp_file.name, 'rb') as f:
                response = HttpResponse(f.read(), content_type='application/zip')
                response['Content-Disposition'] = f'attachment; filename="gene_analysis_results_{task_id}.zip"'
                
            os.unlink(tmp_file.name)
            
            return response
            
    except GeneAnalysisTask.DoesNotExist:
        return JsonResponse({'error': '任务不存在'}, status=404)

# ===================== 辅助函数 =====================
def validate_sequence(sequence):
    valid_bases = {'A', 'T', 'C', 'G', 'a', 't', 'c', 'g'}
    return all(base in valid_bases for base in sequence)

def perform_analysis(sequence):
    length = len(sequence)
    
    counts = {
        'A': sequence.upper().count('A'),
        'T': sequence.upper().count('T'),
        'C': sequence.upper().count('C'),
        'G': sequence.upper().count('G')
    }
    
    gc_content = (counts['G'] + counts['C']) / length * 100 if length > 0 else 0
    
    orfs = find_orfs(sequence)
    
    results = {
        'sequence_length': length,
        'base_composition': counts,
        'gc_content': round(gc_content, 2),
        'orfs': orfs,
        'timestamp': timezone.now().strftime('%Y-%m-%d %H:%M:%S')
    }
    
    return results

def find_orfs(sequence):
    orfs = []
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    
    start_positions = [i for i in range(len(sequence)-2) if sequence[i:i+3] == start_codon]
    
    for start in start_positions:
        for i in range(start, len(sequence)-2, 3):
            codon = sequence[i:i+3]
            if codon in stop_codons:
                orfs.append({
                    'start': start,
                    'end': i+2,
                    'length': i+3 - start,
                    'sequence': sequence[start:i+3]
                })
                break
    
    return orfs

def save_analysis(user, sequence, results):
    try:
        history = AnalysisHistory.objects.create(
            user=user,
            sequence=sequence[:5000],
            results=json.dumps(results),
            created_at=timezone.now()
        )
        return history
    except Exception as e:
        logger.error(f"保存分析历史失败: {str(e)}")
        return None

def biochat_view(request):
    return render(request, "hormone_app/ai.html")
    
@csrf_exempt  # <-- 添加这一行
def get_biogpt_response(request):
    if request.method != "POST":
        return JsonResponse(
            {"status": "error", "message": "Only POST requests are supported"}, status=405
        )

    try:
        # 解析 JSON 请求体
        try:
            request_data = json.loads(request.body)
        except json.JSONDecodeError:
            return JsonResponse(
                {"status": "error", "message": "Invalid JSON format in request body"}, status=400
            )

        input_text = request_data.get("input_text", "").strip()
        if not input_text:
            return JsonResponse(
                {"status": "error", "message": "Please enter question content"}, status=400
            )

        # ========== 核心：混合检索 + 智能响应 ==========
        retrieval_result = hybrid_advanced_retrieval(input_text)
        prompt_or_result = build_biogpt_prompt(input_text, retrieval_result)

        # --- 情况1: 短词查询 → 直接返回结构化数据 ---
        if prompt_or_result["type"] == "structured":
            return JsonResponse({
                "status": "success",
                "response": prompt_or_result["content"],
                "source": "database"
            })

        # --- 情况2: 长句查询 → 调用 BioGPT 生成 ---
        elif prompt_or_result["type"] == "generated":
            final_prompt = prompt_or_result["prompt"]
            
            # 构建发送给 BioGPT API 的请求数据
            api_request_data = {
                "full_prompt": final_prompt,
                "max_new_tokens": 200,
                "max_input_length": 2000,
                "truncation": True,
            }

            try:
                api_response = requests.post(
                    BIOGPT_API_URL,
                    json=api_request_data,
                    timeout=240,
                    headers={"Content-Type": "application/json"}
                )
                api_response.raise_for_status()
                result = api_response.json()

                if result.get("status") == "success" and "response" in result:
                    raw_response = result["response"]
                    
                    # 清理模型可能的前缀（如重复 Prompt）
                    if final_prompt in raw_response:
                        generated_text = raw_response.replace(final_prompt, "").strip()
                    else:
                        generated_text = raw_response.strip()
                    
                    # 兜底：如果模型胡说，强制使用检索结果（可选）
                    bad_phrases = ["it is important", "you have an understanding", "in order to understand"]
                    if any(phrase in generated_text.lower() for phrase in bad_phrases):
                        # 从语义检索中提取激素名（通用版）
                        hormone_names = set()
                        for res in retrieval_result.get("semantic_content", []):
                            match = re.search(r'Hormone ([^i][^s][^ ]+.*?) is associated with', res['text'])
                            if match:
                                hormone_names.add(match.group(1).strip())
                        if hormone_names:
                            generated_text = f"The animal hormones related to the query are: {', '.join(sorted(hormone_names))}."

                    return JsonResponse({
                        "status": "success",
                        "response": generated_text if generated_text else "No response generated.",
                        "source": "database" if retrieval_result["has_db_data"] else "model"
                    })
                else:
                    return JsonResponse({
                        "status": "error",
                        "message": "BioGPT returned invalid response"
                    }, status=502)

            except requests.RequestException as req_err:
                return JsonResponse({
                    "status": "error",
                    "message": f"BioGPT API request failed: {str(req_err)[:150]}"
                }, status=502)

        else:
            return JsonResponse({
                "status": "error",
                "message": "Unknown response type from build_biogpt_prompt"
            }, status=500)

    except Exception as e:
        logger.error(f"get_biogpt_response error: {str(e)}", exc_info=True)
        return JsonResponse(
            {"status": "error", "message": f"Internal server error: {str(e)[:150]}"},
            status=500
        )
@login_required
def analysis(request):
    if request.method == 'POST':
        form = AnalysisJobForm(request.POST)
        if form.is_valid():
            job = form.save(commit=False)
            job.user = request.user
            job.save()
            return redirect('analysis_detail', job_id=job.id)
    else:
        form = AnalysisJobForm()
    
    jobs = AnalysisJob.objects.filter(user=request.user).order_by('-created_at')[:5]
    
    return render(request, 'hormone_app/docs.html', {
        'form': form,
        'jobs': jobs
    })

@login_required
def analysis_detail(request, job_id):
    job = get_object_or_404(AnalysisJob, id=job_id, user=request.user)
    
    uploaded_files = job.uploaded_files.all()
    file_types = {f.file_type: f for f in uploaded_files}
    
    return render(request, 'hormone_app/analysis_detail.html', {
        'job': job,
        'file_types': file_types,
        'upload_form': FileUploadForm()
    })

@login_required
def upload_file(request, job_id):
    job = get_object_or_404(AnalysisJob, id=job_id, user=request.user)
    
    if request.method == 'POST' and request.FILES.get('file'):
        form = FileUploadForm(request.POST, request.FILES)
        if form.is_valid():
            file_type = request.POST.get('file_type')
            uploaded_file = request.FILES['file']
            
            file = form.save(commit=False)
            file.job = job
            file.file_type = file_type
            file.filename = uploaded_file.name
            file.file_size = uploaded_file.size
            file.save()
            
            return JsonResponse({
                'status': 'success',
                'message': f'File "{uploaded_file.name}" uploaded successfully',
                'file_type': file_type,
                'filename': uploaded_file.name
            })
    
    return JsonResponse({'status': 'error', 'message': 'File upload failed'}, status=400)

@login_required
def save_parameters(request, job_id):
    job = get_object_or_404(AnalysisJob, id=job_id, user=request.user)
    
    if request.method == 'POST':
        try:
            data = json.loads(request.body)
            parameter_group = data.get('group')
            parameters = data.get('parameters', {})
            
            for name, value in parameters.items():
                if isinstance(value, list):
                    value = ','.join(value)
                
                param, created = AnalysisParameter.objects.get_or_create(
                    job=job,
                    group=parameter_group,
                    name=name,
                    defaults={'value': str(value)}
                )
                
                if not created:
                    param.value = str(value)
                    param.save()
            
            return JsonResponse({'status': 'success', 'message': 'Parameters saved'})
        except Exception as e:
            return JsonResponse({'status': 'error', 'message': str(e)}, status=400)
    
    return JsonResponse({'status': 'error', 'message': 'Invalid request method'}, status=405)

@login_required
def run_analysis(request, job_id):
    job = get_object_or_404(AnalysisJob, id=job_id, user=request.user)
    
    required_files = job.uploaded_files.filter(file_type__in=['methyl_data', 'metadata'])
    if required_files.count() < 2:
        return JsonResponse({
            'status': 'error',
            'message': 'Please upload both methylation data and metadata files'
        }, status=400)
    
    job.status = 'running'
    job.started_at = timezone.now()
    job.progress = 0
    job.save()
    
    run_analysis_pipeline.delay(str(job.id))
    
    return JsonResponse({
        'status': 'success',
        'message': 'Analysis started successfully'
    })

@login_required
def check_progress(request, job_id):
    job = get_object_or_404(AnalysisJob, id=job_id, user=request.user)
    
    return JsonResponse({
        'status': job.status,
        'progress': job.progress,
        'log': job.log
    })

@login_required
def get_results(request, job_id):
    job = get_object_or_404(AnalysisJob, id=job_id, user=request.user)
    
    if job.status != 'completed':
        return JsonResponse({
            'status': 'error',
            'message': 'Analysis not completed yet'
        }, status=400)
    
    results = job.results.all()
    result_data = []
    
    for result in results:
        result_data.append({
            'id': str(result.id),
            'type': result.result_type,
            'type_display': result.get_result_type_display(),
            'filename': result.filename,
            'description': result.description,
            'created_at': result.created_at.strftime('%Y-%m-%d %H:%M:%S')
        })
    
    stats = {
        'total_cpgs': '485,577',
        'total_dmps': '1,243',
        'total_dmrs': '87',
        'associated_genes': '342'
    }
    
    return JsonResponse({
        'status': 'success',
        'results': result_data,
        'statistics': stats
    })

@login_required
def download_result(request, result_id):
    result = get_object_or_404(AnalysisResult, id=result_id, job__user=request.user)
    
    with open(result.file.path, 'rb') as f:
        response = HttpResponse(f.read(), content_type='application/octet-stream')
        response['Content-Disposition'] = f'attachment; filename="{result.filename}"'
        return response

@login_required
def get_visualization_data(request, job_id, plot_type):
    job = get_object_or_404(AnalysisJob, id=job_id, user=request.user)
    
    if job.status != 'completed':
        return JsonResponse({
            'status': 'error',
            'message': 'Analysis not completed yet'
        }, status=400)
    
    sample_data = {
        'methylation-distribution': {
            'x': [i/100 for i in range(101)],
            'control': [0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1],
            'treatment': [0.2, 0.3, 0.5, 0.7, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05]
        },
        'pca-analysis': {
            'pc1': [1.2, 1.5, 1.3, 1.4, -1.1, -1.3, -1.2, -1.4],
            'pc2': [0.8, 0.7, 0.9, 0.8, -0.7, -0.8, -0.9, -0.7],
            'groups': ['control', 'control', 'control', 'control', 'treatment', 'treatment', 'treatment', 'treatment']
        }
    }
    
    data = sample_data.get(plot_type, {})
    
    return JsonResponse({
        'status': 'success',
        'data': data
    })

def is_superuser(user):
    return user.is_superuser

@login_required
@user_passes_test(is_superuser)
def get_model_data(request):
    model_name = request.GET.get('model')
    if not model_name:
        return JsonResponse({'success': False, 'error': '未指定模型名称'})
    
    try:
        model = apps.get_model('hormone_app', model_name)
        records = model.objects.all().values()
        return JsonResponse({
            'success': True,
            'data': list(records)
        })
    except Exception as e:
        return JsonResponse({'success': False, 'error': str(e)})
        
        
        
# --------------------------
# 1. Smith-Waterman 局部比对算法（DNA序列优化）
# --------------------------
class SmithWaterman:
    """Smith-Waterman局部序列比对算法（适用于DNA序列）"""
    def __init__(self, match_score=2, mismatch_penalty=-1, gap_penalty=-1):
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty

    def align(self, seq1, seq2, max_length=1000):
        """
        执行Smith-Waterman比对
        :param seq1: 输入序列
        :param seq2: 数据库序列
        :param max_length: 长序列截断阈值（优化速度）
        :return: 比对结果 + 相似度
        """
        # 长序列截断优化（仅比对前max_length个碱基）
        seq1_trunc = seq1[:max_length].upper() if len(seq1) > max_length else seq1.upper()
        seq2_trunc = seq2[:max_length].upper() if len(seq2) > max_length else seq2.upper()

        # 初始化得分矩阵
        rows = len(seq1_trunc) + 1
        cols = len(seq2_trunc) + 1
        score_matrix = [[0 for _ in range(cols)] for _ in range(rows)]
        traceback_matrix = [[0 for _ in range(cols)] for _ in range(rows)]

        # 填充得分矩阵
        max_score = 0
        max_pos = (0, 0)

        for i in range(1, rows):
            for j in range(1, cols):
                # 计算匹配/错配得分
                if seq1_trunc[i-1] == seq2_trunc[j-1]:
                    match = score_matrix[i-1][j-1] + self.match_score
                else:
                    match = score_matrix[i-1][j-1] + self.mismatch_penalty

                # 计算缺口得分
                delete = score_matrix[i-1][j] + self.gap_penalty
                insert = score_matrix[i][j-1] + self.gap_penalty

                # 取最大值（最小为0）
                score = max(match, delete, insert, 0)
                score_matrix[i][j] = score

                # 记录回溯方向
                if score == 0:
                    traceback_matrix[i][j] = 0  # 停止
                elif score == match:
                    traceback_matrix[i][j] = 1  # 对角线
                elif score == delete:
                    traceback_matrix[i][j] = 2  # 上
                else:
                    traceback_matrix[i][j] = 3  # 左

                # 更新最大得分位置
                if score > max_score:
                    max_score = score
                    max_pos = (i, j)

        # 回溯获取比对结果
        align1 = []
        align2 = []
        match_line = []
        i, j = max_pos

        while traceback_matrix[i][j] != 0:
            if traceback_matrix[i][j] == 1:
                # 匹配/错配
                align1.append(seq1_trunc[i-1])
                align2.append(seq2_trunc[j-1])
                if seq1_trunc[i-1] == seq2_trunc[j-1]:
                    match_line.append('|')  # 匹配
                else:
                    match_line.append('*')  # 错配
                i -= 1
                j -= 1
            elif traceback_matrix[i][j] == 2:
                # 序列1缺口
                align1.append(seq1_trunc[i-1])
                align2.append('-')
                match_line.append(' ')
                i -= 1
            elif traceback_matrix[i][j] == 3:
                # 序列2缺口
                align1.append('-')
                align2.append(seq2_trunc[j-1])
                match_line.append(' ')
                j -= 1

        # 反转得到正确顺序
        align1 = ''.join(reversed(align1))
        align2 = ''.join(reversed(align2))
        match_line = ''.join(reversed(match_line))

        # 计算相似度（基于比对区域）
        match_count = match_line.count('|')
        total_aligned = len(match_line)
        similarity = (match_count / total_aligned * 100) if total_aligned > 0 else 0

        return {
            'similarity': round(similarity, 2),
            'alignment': {
                'input_seq': align1,
                'db_seq': align2,
                'match_line': match_line
            },
            'max_score': max_score,
            'truncated': len(seq1) > max_length or len(seq2) > max_length
        }

# --------------------------
# 2. 基因预测主视图
# --------------------------


@csrf_protect
def gene_prediction(request):
    if request.method == 'GET':
        return render(request, 'prediction.html')
    
    elif request.method == 'POST':
        try:
            # 解析前端数据
            data = json.loads(request.body)
            prediction_type = data.get('prediction_type', 'hormone')
            sequences = data.get('sequences', [])
            species = data.get('species', '')
            threshold = int(data.get('similarity_threshold', 50))

            # 验证输入
            if not sequences:
                return JsonResponse({
                    'success': False,
                    'total_matched': 0,
                    'results': {},
                    'error': 'No sequences provided'
                }, status=400)

            # 初始化Smith-Waterman比对器
            sw = SmithWaterman(match_score=2, mismatch_penalty=-1, gap_penalty=-1)

            # 查询所有数据库序列
            all_records = HormoneRelatedGene.objects.all().values(
                'id', 'hormone_name', 'related_genes', 'pmid', 'gene_sequence', 'related_diseases'
            )

            # 批量比对所有输入序列
            results = {}
            total_matched = 0

            for input_seq in sequences:
                seq_header = input_seq.get('header', 'Unknown')
                seq_str = input_seq.get('sequence', '')
                
                if not seq_str:
                    results[seq_header] = []
                    continue

                # 比对当前输入序列与所有数据库序列
                seq_matched = []
                for db_record in all_records:
                    db_seq = db_record.get('gene_sequence', '')
                    if not db_seq:
                        continue

                    # 执行Smith-Waterman比对
                    align_result = sw.align(seq_str, db_seq, max_length=1000)  # 长序列截断为1000bp
                    similarity = align_result['similarity']

                    # 筛选阈值以上的结果
                    if similarity >= threshold:
                        db_record['similarity'] = similarity
                        db_record['alignment'] = align_result['alignment']
                        db_record['truncated'] = align_result['truncated']
                        seq_matched.append(db_record)
                        total_matched += 1

                # 按相似度降序排序
                seq_matched.sort(key=lambda x: x['similarity'], reverse=True)
                results[seq_header] = seq_matched

            return JsonResponse({
                'success': True,
                'total_matched': total_matched,
                'results': results,
                'species': species,
                'prediction_type': prediction_type,
                'threshold': threshold
            })
        
        except Exception as e:
            logger.error(f"基因预测失败: {str(e)}", exc_info=True)  # 增加日志记录
            return JsonResponse({
                'success': False,
                'total_matched': 0,
                'results': {},
                'error': str(e)
            }, status=500)

# --------------------------
# 3. PDF报告下载视图（精简表格版）
# --------------------------
@csrf_protect
def download_pdf(request):
    if request.method == 'POST':
        try:
            # 获取表单数据并验证
            prediction_type = request.POST.get('prediction_type', 'hormone')
            species = request.POST.get('species', 'Unknown').strip() or 'Unknown'  # 处理空值
            threshold = request.POST.get('similarity_threshold', 50)
            
            # 解析JSON数据
            try:
                matched_results = json.loads(request.POST.get('matched_results', '{}'))
            except json.JSONDecodeError as e:
                logger.error(f"匹配结果JSON解析失败: {str(e)}")
                return JsonResponse({
                    'success': False,
                    'error': 'Invalid data format for results'
                }, status=400)

            # 创建PDF响应
            response = HttpResponse(content_type='application/pdf')
            filename = f'gene_prediction_report_{datetime.now().strftime("%Y%m%d_%H%M%S")}.pdf'
            response['Content-Disposition'] = (
                f'attachment; filename="{quote(filename)}"; '
                f'filename*=UTF-8\'\'{quote(filename)}'
            )

            # 初始化PDF文档（调整列宽适配删除后的表格）
            doc = SimpleDocTemplate(
                response,
                pagesize=A4,
                rightMargin=20,
                leftMargin=20,
                topMargin=30,
                bottomMargin=20
            )
            styles = getSampleStyleSheet()
            elements = []

            # --------------------------
            # 标题和实验室信息（不变）
            # --------------------------
            title_style = ParagraphStyle(
                'CustomTitle',
                parent=styles['Heading1'],
                fontSize=18,
                spaceAfter=20,
                alignment=1,
                textColor=colors.darkblue
            )
            elements.append(Paragraph('Gene Prediction Report (Smith-Waterman Alignment)', title_style))
            elements.append(Spacer(1, 12))

            lab_style = ParagraphStyle(
                'LabInfo',
                parent=styles['Normal'],
                fontSize=11,
                spaceAfter=15,
                alignment=1
            )
            elements.append(Paragraph('Inner Mongolia University of Science and Technology', lab_style))
            elements.append(Paragraph('Animal Hormone Research Lab', lab_style))
            elements.append(Paragraph(f'Contact: contact@hormone-db.com | Phone: +86 15393345650', lab_style))
            elements.append(Spacer(1, 20))

            # --------------------------
            # 基本信息（不变）
            # --------------------------
            info_style = ParagraphStyle(
                'Info',
                parent=styles['Normal'],
                fontSize=12,
                spaceAfter=10
            )
            elements.append(Paragraph(f'Report Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}', info_style))
            elements.append(Paragraph(f'Prediction Type: {prediction_type.capitalize()} Prediction', info_style))
            elements.append(Paragraph(f'Species: {species}', info_style))
            elements.append(Paragraph(f'Similarity Threshold: ≥ {threshold}%', info_style))
            elements.append(Spacer(1, 20))

            # --------------------------
            # 结果表格（删除Alignment Preview列）
            # --------------------------
            if matched_results:
                for seq_header, results in matched_results.items():
                    seq_title = Paragraph(f'Sequence: {seq_header}', styles['Heading2'])
                    elements.append(seq_title)
                    elements.append(Spacer(1, 10))

                    if results:
                        # 表格数据（移除Alignment Preview列）
                        table_data = [
                            ['Similarity (%)', 'Hormone/Disease', 'Related Info', 'PubMed ID']  # 列头精简
                        ]

                        for result in results:
                            if prediction_type == 'hormone':
                                main_content = result.get('hormone_name', 'N/A')
                                sub_content = result.get('related_genes', 'N/A')
                            else:
                                main_content = 'Related Diseases'
                                sub_content = result.get('related_diseases', 'N/A')

                            # 表格行数据（不再包含Alignment Preview）
                            table_data.append([
                                f"{result.get('similarity', 0)}%",
                                main_content,
                                sub_content,
                                str(result.get('pmid', 'N/A'))
                            ])

                        # 调整表格列宽（适配4列，更宽松）
                        table = Table(
                            table_data,
                            colWidths=[1.2*inch, 2.2*inch, 2.8*inch, 1.5*inch]  # 4列宽度分配
                        )
                        table_style = TableStyle([
                            ('BACKGROUND', (0, 0), (-1, 0), colors.lightblue),
                            ('TEXTCOLOR', (0, 0), (-1, 0), colors.black),
                            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                            ('FONTSIZE', (0, 0), (-1, 0), 10),
                            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                            ('BACKGROUND', (0, 1), (-1, -1), colors.white),
                            ('GRID', (0, 0), (-1, -1), 1, colors.black),
                            ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
                            ('FONTSIZE', (0, 1), (-1, -1), 9),
                            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                            ('WORDWRAP', (0, 0), (-1, -1), 'CJK'),
                        ])
                        table.setStyle(table_style)
                        elements.append(table)
                        elements.append(Spacer(1, 20))

                        # 详细序列比对（保留，作为完整展示）
                        align_style = ParagraphStyle(
                            'Alignment',
                            parent=styles['Code'],
                            fontSize=9,
                            spaceAfter=15,
                            fontName='Courier'
                        )
                        elements.append(Paragraph('Detailed Sequence Alignment:', styles['Heading3']))
                        for result in results:
                            hormone_name = result.get('hormone_name', 'Unknown')
                            similarity = result.get('similarity', 0)
                            elements.append(Paragraph(f"> {hormone_name} (Similarity: {similarity}%)", styles['Normal']))
                            alignment = result.get('alignment', {})
                            elements.append(Paragraph(alignment.get('input_seq', ''), align_style))
                            elements.append(Paragraph(alignment.get('match_line', ''), align_style))
                            elements.append(Paragraph(alignment.get('db_seq', ''), align_style))
                            elements.append(Spacer(1, 10))
                    else:
                        elements.append(Paragraph(f'No matches for this sequence (similarity < {threshold}%)', info_style))
                        elements.append(Spacer(1, 10))
            else:
                no_result_style = ParagraphStyle(
                    'NoResult',
                    parent=styles['Normal'],
                    fontSize=12,
                    spaceAfter=10,
                    textColor=colors.red,
                    alignment=1
                )
                elements.append(Paragraph(f'No matching records found (similarity < {threshold}%).', no_result_style))

            # --------------------------
            # 页脚（不变）
            # --------------------------
            footer_style = ParagraphStyle(
                'Footer',
                parent=styles['Normal'],
                fontSize=10,
                spaceAfter=10,
                alignment=1,
                textColor=colors.gray
            )
            elements.append(Spacer(1, 30))
            elements.append(Paragraph('--- End of Report ---', footer_style))
            elements.append(Paragraph('Generated by Animal Hormone Regulation Database', footer_style))

            # 构建PDF文档
            doc.build(elements)
            return response
        
        except Exception as e:
            logger.error(f"PDF生成失败: {str(e)}", exc_info=True)
            return JsonResponse({
                'success': False,
                'error': f'PDF generation failed: {str(e)}'
            }, status=500)
    
    return JsonResponse({
        'success': False,
        'error': 'Method not allowed'
    }, status=405)



from django.db import connection

from django.db import connection as db_connection

def api_hormone_data(request):
    """API: return paginated hormone_data with search support"""
    try:
        page = int(request.GET.get('page', 1))
        page_size = int(request.GET.get('page_size', 20))
        search = request.GET.get('search', '').strip()
        organism_filter = request.GET.get('organism', '').strip()
        
        with db_connection.cursor() as cursor:
            where_clauses = []
            params = []
            
            if search:
                where_clauses.append("(hormone_name LIKE %s OR related_genes LIKE %s OR related_diseases LIKE %s)")
                params.extend([f'%{search}%'] * 3)
            
            if organism_filter:
                where_clauses.append('organism = %s')
                params.append(organism_filter)
            
            where_sql = ' AND '.join(where_clauses) if where_clauses else '1=1'
            
            cursor.execute(f'SELECT COUNT(*) FROM hormone_data WHERE {where_sql}', params)
            total = cursor.fetchone()[0]
            
            offset = (page - 1) * page_size
            sql = f"SELECT id, hormone_name, organism, related_genes, regulation_type, pmid, related_diseases, DO_ID, doid_standardized_name, gene_description, chromosome FROM hormone_data WHERE {where_sql} ORDER BY id LIMIT %s OFFSET %s"
            cursor.execute(sql,
                params + [page_size, offset])
            
            rows = cursor.fetchall()
            data = []
            for row in rows:
                data.append({
                    'id': row[0], 'hormone_name': row[1] or '', 'organism': row[2] or '',
                    'related_genes': row[3] or '', 'regulation_type': row[4] or '',
                    'pmid': str(row[5]) if row[5] else '',
                    'related_diseases': (row[6] or '')[:200],
                    'DO_ID': row[7] or '', 'doid_standardized_name': row[8] or '',
                    'gene_description': (row[9] or '')[:100], 'chromosome': row[10] or '',
                })
            
            cursor.execute('SELECT DISTINCT organism FROM hormone_data WHERE organism IS NOT NULL ORDER BY organism')
            organisms = [r[0] for r in cursor.fetchall()]
            
            return JsonResponse({
                'data': data, 'total': total, 'page': page, 'page_size': page_size,
                'total_pages': (total + page_size - 1) // page_size, 'organisms': organisms
            })
    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)

def api_hormone_stats(request):
    """API: return statistics for charts"""
    try:
        with db_connection.cursor() as cursor:
            cursor.execute('SELECT organism, COUNT(*) as cnt FROM hormone_data GROUP BY organism ORDER BY cnt DESC LIMIT 20')
            species_data = [{'name': r[0], 'value': r[1]} for r in cursor.fetchall()]
            
            cursor.execute('SELECT regulation_type, COUNT(*) as cnt FROM hormone_data GROUP BY regulation_type')
            regulation_data = [{'name': r[0] or 'Unknown', 'value': r[1]} for r in cursor.fetchall()]
            
            cursor.execute('SELECT related_diseases, COUNT(*) as cnt FROM hormone_data WHERE related_diseases IS NOT NULL AND related_diseases != "" GROUP BY related_diseases ORDER BY cnt DESC LIMIT 15')
            disease_data = [{'name': (r[0] or '')[:50], 'value': r[1]} for r in cursor.fetchall()]
            
            cursor.execute('SELECT COUNT(*), COUNT(DISTINCT hormone_name), COUNT(DISTINCT related_genes), COUNT(DISTINCT organism) FROM hormone_data')
            row = cursor.fetchone()
            summary = {'total_records': row[0], 'unique_hormones': row[1], 'unique_genes': row[2], 'unique_species': row[3]}
            
            return JsonResponse({'species': species_data, 'regulation': regulation_data, 'diseases': disease_data, 'summary': summary})
    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)

# ==================== Chart Data API ====================
@csrf_exempt
def chart_data(request):
    """返回主页图表所需的真实数据"""
    from hormone_app.models import HormoneRelatedGene
    from django.db.models import Count
    
    # 1. Top 10物种分布
    species = list(HormoneRelatedGene.objects.exclude(organism='').values('organism').annotate(count=Count('id')).order_by('-count')[:10])
    
    # 2. 调控类型分布
    reg_types = list(HormoneRelatedGene.objects.exclude(regulation_type='').values('regulation_type').annotate(count=Count('id')).order_by('-count'))
    
    # 3. Top 8疾病
    diseases = list(HormoneRelatedGene.objects.exclude(related_diseases='').exclude(related_diseases='Not available').values('related_diseases').annotate(count=Count('id')).order_by('-count')[:8])
    
    # 4. Top 10激素
    hormones = list(HormoneRelatedGene.objects.values('hormone_name').annotate(count=Count('id')).order_by('-count')[:10])
    
    # 5. 物种-疾病热力图数据 (Top 5物种 × Top 6疾病)
    top5_species = [s['organism'] for s in species[:5]]
    top6_diseases = [d['related_diseases'] for d in diseases[:6]]
    heatmap = []
    for di, disease in enumerate(top6_diseases):
        for si, sp in enumerate(top5_species):
            cnt = HormoneRelatedGene.objects.filter(organism=sp, related_diseases=disease).count()
            if cnt > 0:
                heatmap.append([si, di, cnt])
    
    return JsonResponse({
        'species': species,
        'reg_types': reg_types,
        'diseases': diseases,
        'hormones': hormones,
        'heatmap': {
            'species': top5_species,
            'diseases': top6_diseases,
            'data': heatmap
        }
    })

# ---- PanHorm v2.6.3 safe sequence release-index API (read-only, added 2026-05-19) ----
def api_v263_sequence_index(request):
    """Read-only API for the v2.6.3 safe sequence release-index table.

    Query params:
      q: keyword matched against hormone/gene/disease/DO fields
      gene: exact/partial standard gene symbol or NCBI gene id
      hormone: partial hormone name
      species: partial species
      status: v2_6_3_sequence_merge_status
      page, page_size
    """
    from django.db import connection
    from django.views.decorators.http import require_GET
    # Kept inside function to avoid import side effects in this large legacy views.py.
    try:
        page = max(1, int(request.GET.get("page", 1)))
        page_size = min(100, max(1, int(request.GET.get("page_size", 20))))
    except ValueError:
        return JsonResponse({"ok": False, "error": "invalid page/page_size"}, status=400)

    allowed_status = {"", "merged_pass", "not_merged_warn", "not_merged_fail", "no_sequence_candidate"}
    status_filter = request.GET.get("status", "").strip()
    if status_filter not in allowed_status:
        return JsonResponse({"ok": False, "error": "invalid status"}, status=400)

    secondary_status_filter = request.GET.get("secondary_status", "").strip()
    allowed_secondary = {"", "PASS", "REVIEW", "FAIL", "NOT_MERGED"}
    if secondary_status_filter not in allowed_secondary:
        return JsonResponse({"ok": False, "error": "invalid secondary_status"}, status=400)

    q = request.GET.get("q", "").strip()
    gene = request.GET.get("gene", "").strip()
    hormone = request.GET.get("hormone", "").strip()
    species = request.GET.get("species", "").strip()

    where = []
    params = []
    if status_filter:
        where.append("v2_6_3_sequence_merge_status=%s")
        params.append(status_filter)
    if secondary_status_filter:
        where.append("secondary_qc_status=%s")
        params.append(secondary_status_filter)
    if gene:
        where.append("(standard_gene_symbol LIKE %s OR original_gene_symbol LIKE %s OR ncbi_gene_id LIKE %s)")
        like = f"%{gene}%"
        params.extend([like, like, like])
    if hormone:
        where.append("hormone_name LIKE %s")
        params.append(f"%{hormone}%")
    if species:
        where.append("standard_species LIKE %s")
        params.append(f"%{species}%")
    if q:
        like = f"%{q}%"
        where.append("(hormone_name LIKE %s OR standard_gene_symbol LIKE %s OR disease_name LIKE %s OR do_term LIKE %s OR doid LIKE %s)")
        params.extend([like, like, like, like, like])

    where_sql = " WHERE " + " AND ".join(where) if where else ""
    offset = (page - 1) * page_size
    fields = [
        "source_csv_row", "hormone_name", "hormone_accession", "standard_gene_symbol", "original_gene_symbol",
        "ncbi_gene_id", "standard_species", "tax_id", "disease_name", "doid", "do_term",
        "gene_disease_evidence_level", "gene_disease_source", "hg_pmid", "gd_pmid",
        "protein_accession", "mRNA_accession", "sequence_match_level", "sequence_source",
        "v2_6_3_sequence_merge_status", "v2_6_3_sequence_qc_severity", "v2_6_3_sequence_qc_flags",
        "v2_6_3_ncbi_mRNA_accession", "v2_6_3_ncbi_protein_accession", "v2_6_3_ncbi_gene_tax_id",
        "v2_6_3_ncbi_locus", "v2_6_3_sequence_source", "v2_6_3_sequence_source_gene_id",
        "secondary_qc_status", "secondary_qc_flags", "secondary_recommended_action"
    ]
    try:
        with connection.cursor() as cursor:
            cursor.execute(f"SELECT COUNT(*) FROM panhorm_v263_safe_sequence_index{where_sql}", params)
            total = cursor.fetchone()[0]
            cursor.execute(
                "SELECT " + ",".join(f"`{x}`" for x in fields) +
                f" FROM panhorm_v263_safe_sequence_index{where_sql} ORDER BY id LIMIT %s OFFSET %s",
                params + [page_size, offset]
            )
            rows = [dict(zip(fields, row)) for row in cursor.fetchall()]
    except Exception as exc:
        return JsonResponse({"ok": False, "error": str(exc)[:300]}, status=500)
    return JsonResponse({
        "ok": True,
        "version": "v2.6.3_safe_sequences_index",
        "note": "Read-only release index for QC-passed/candidate sequence annotations. Full long sequences are stored in the controlled CSV sidecar under /root/panhorm_v263_import_20260519 and are not exposed by this API.",
        "page": page,
        "page_size": page_size,
        "total": total,
        "results": rows,
    }, json_dumps_params={"ensure_ascii": False})


def api_v263_sequence_stats(request):
    """Read-only summary stats for PanHorm v2.6.3 safe sequence release index."""
    from django.db import connection
    try:
        with connection.cursor() as cursor:
            cursor.execute("SELECT COUNT(*) FROM panhorm_v263_safe_sequence_index")
            total = cursor.fetchone()[0]
            cursor.execute("SELECT v2_6_3_sequence_merge_status, COUNT(*) FROM panhorm_v263_safe_sequence_index GROUP BY v2_6_3_sequence_merge_status")
            status_counts = dict(cursor.fetchall())
            cursor.execute("SELECT standard_species, COUNT(*) c FROM panhorm_v263_safe_sequence_index GROUP BY standard_species ORDER BY c DESC LIMIT 20")
            species_counts = [{"species": r[0], "count": r[1]} for r in cursor.fetchall()]
            cursor.execute("SELECT secondary_qc_status, COUNT(*) FROM panhorm_v263_safe_sequence_index GROUP BY secondary_qc_status")
            secondary_status_counts = dict(cursor.fetchall())
    except Exception as exc:
        return JsonResponse({"ok": False, "error": str(exc)[:300]}, status=500)
    return JsonResponse({
        "ok": True,
        "version": "v2.6.3_safe_sequences_index",
        "total": total,
        "status_counts": status_counts,
        "secondary_status_counts": secondary_status_counts,
        "top_species": species_counts,
    }, json_dumps_params={"ensure_ascii": False})
# ---- End PanHorm v2.6.3 safe sequence release-index API ----
