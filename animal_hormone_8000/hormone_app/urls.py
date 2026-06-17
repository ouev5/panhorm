from django.urls import path
from django.views.decorators.csrf import csrf_exempt
from django.contrib.auth.views import LogoutView
from . import views
from .analysishou import see_view

from hormone_app.analysishou import gnn_views, phylogeny_views 
# 正确导入分析视图 - 从当前目录下的analysishou子目录导入
from .analysishou.sc_analysis_view import analysis_page, create_task, upload_file, run_analysis, get_task_status, get_task_results, download_result, download_all_results
# --- 修改导入部分：增加新视图函数 ---
from .analysishou.sc_analysis_view import get_visualization_plot, download_result_by_type
from .analysishou.sc_analysis_view import gene_expression_analysis  # 新增基因表达分析视图导入

# --- 新增导入：导入 BLAST 视图函数 ---
from .analysishou.blast_view import blast_create_task, blast_run_analysis, blast_get_task_status, blast_get_task_results, blast_download_result_csv, blast_download_full_report_xml, blast_download_fasta_hits, blast_download_all_results, blast_download_evalue_plot, blast_download_identity_plot, blast_download_all_visualizations,blast_get_alignment_details

# --- 新增导入：导入 DNA甲基化分析视图函数 ---
from .analysishou.DNA_Methylation_view import dna_methylation_analysis_page, dna_create_task, dna_upload_file, dna_run_analysis, dna_get_task_status, dna_get_task_results, dna_download_result, dna_download_all_results, dna_get_visualization_plot, dna_download_result_by_type, test_dna_url

# --- 新增导入：导入转录组分析视图函数 ---
from .analysishou.transcriptome_view import (
    transcriptome_analysis_page,
    create_transcriptome_task,
    upload_transcriptome_file,
    run_transcriptome_analysis,
    get_transcriptome_task_status,
    get_transcriptome_task_results,
    get_transcriptome_visualization_plot,
    download_transcriptome_result_by_type,
    download_all_transcriptome_results,
    list_user_transcriptome_tasks,
    cancel_transcriptome_task
)

# 在导入部分添加KEGG视图
from .analysishou.Kegg_view import (
    kegg_analysis_page,
    create_kegg_task,
    upload_kegg_file,
    run_kegg_analysis,
    get_kegg_task_status,
    get_kegg_task_results,
    # 移除 download_kegg_result，使用 download_kegg_result_by_type
    download_all_kegg_results,
    get_kegg_visualization_plot,
    download_kegg_result_by_type,
    list_user_kegg_tasks,
    cancel_kegg_task,
    validate_kegg_gene_list,
    get_kegg_previous_results,
    save_kegg_parameters,
    download_kegg_analysis_report,
    preview_kegg_mapping,
    get_kegg_database_stats,
    get_kegg_pathway_detail
)

from .views import UnifiedDeepSeekView  # 导入视图类


app_name = 'hormone_app'

urlpatterns = [
    # 认证相关
    path('login_view/', views.login_view, name='login_view'),
    path('register/', views.register_view, name='register'),
    path('logout/', LogoutView.as_view(next_page='/'), name='logout'),

    # 基础页面
    path('', views.home, name='home'),
    path('docs/', views.docs, name='docs'),
    path('data/', views.data, name='data'),
    path('hormone_list/', views.list_view, name='hormone_list'),
    path('help/', views.help_view, name='help'),
    path('ai/', views.ai_view, name='ai'),  # 前端页面
    path("api/get-response/", views.get_biogpt_response, name="get_biogpt_response"),  # 后端接口

    # 搜索
    path('search/', views.search, name='search'),
    path('entity/<str:entity_type>/', views.entity_detail, name='entity_detail'),

    # 激素受体列表
    path('receptors/', views.hormone_receptor_list, name='receptor_list'),

    # 基因分析相关
    path('gene_analysis/', views.gene_analysis, name='gene_analysis'),
    path('analysis_history/', views.analysis_history, name='analysis_history'),
    path('analysis_detail/<int:analysis_id>/', views.analysis_detail, name='analysis_detail'),

    # 管理员功能
    path('admin/dashboard/', views.admin_dashboard, name='admin_dashboard'),
    path('admin/add_user/', views.add_user, name='add_user'),
    path('admin/delete_user/<int:user_id>/', views.delete_user, name='delete_user'),
    path('admin/upload_literature/', views.upload_literature, name='upload_literature'),
    path('get-model-data/', views.get_model_data, name='get_model_data'),
    
    # 分析页面
    path('docs/1.html', views.page1, name='page1'),
    path('docs/2.html', views.page2, name='page2'),
    path('docs/3.html', views.page3, name='page3'),
    path('docs/4.html', views.page4, name='page4'),
    path('docs/5.html', views.page5, name='page5'),
    path('docs/6.html', views.page6, name='page6'),

    # 单细胞分析功能相关URL - 修正引用方式
    path('analysis/', analysis_page, name='sc_analysis_page'),
    path('api/tasks/', create_task, name='create_analysis_task'),
    path('api/tasks/<uuid:task_id>/upload/', upload_file, name='upload_analysis_file'),
    path('api/tasks/<uuid:task_id>/run/', run_analysis, name='run_analysis_task'),
    path('api/tasks/<uuid:task_id>/status/', get_task_status, name='get_analysis_status'),
    path('api/tasks/<uuid:task_id>/results/', get_task_results, name='get_analysis_results'),
    path('api/results/<uuid:result_id>/download/', download_result, name='download_analysis_result'),
    path('api/tasks/<uuid:task_id>/download-all/', download_all_results, name='download_all_analysis_results'),

    # --- 新增路由：单细胞分析辅助功能 ---
    path('api/tasks/<uuid:task_id>/plots/<str:plot_type>/', get_visualization_plot, name='get_visualization_plot'),
    path('api/tasks/<uuid:task_id>/download-by-type/', download_result_by_type, name='download_result_by_type'),

    # 保留旧的URL兼容（如果需要）
    path('analysis/<uuid:job_id>/', views.analysis_detail, name='analysis_detail'),
    path('analysis/<uuid:job_id>/save-params/', views.save_parameters, name='save_parameters'),
    path('analysis/<uuid:job_id>/progress/', views.check_progress, name='check_progress'),
    path('analysis/<uuid:job_id>/visualization/<str:plot_type>/', views.get_visualization_data, name='get_visualization_data'),

    # --- 新增路由：BLAST API 路由 ---
    # 任务管理
    path('api/blast-tasks/', blast_create_task, name='blast_create_task'),
    # 分析执行
    path('api/blast-tasks/<uuid:task_id>/run/', blast_run_analysis, name='blast_run_analysis'),
    # 状态查询
    path('api/blast-tasks/<uuid:task_id>/status/', blast_get_task_status, name='blast_get_task_status'),
    # 结果获取
    path('api/blast-tasks/<uuid:task_id>/results/', blast_get_task_results, name='blast_get_task_results'),
    path('api/blast-tasks/<uuid:task_id>/alignment/', blast_get_alignment_details, name='blast_alignment'),
    # 结果下载
    path('api/blast-tasks/<uuid:task_id>/download/blast-results/', blast_download_result_csv, name='blast_download_result_csv'),
    path('api/blast-tasks/<uuid:task_id>/download/full-report/', blast_download_full_report_xml, name='blast_download_full_report_xml'),
    path('api/blast-tasks/<uuid:task_id>/download/fasta-hits/', blast_download_fasta_hits, name='blast_download_fasta_hits'),
    path('api/blast-tasks/<uuid:task_id>/download-all/', blast_download_all_results, name='blast_download_all_results'),
    # 可视化图表下载
    path('api/blast-tasks/<uuid:task_id>/download/evalue-plot/', blast_download_evalue_plot, name='blast_download_evalue_plot'),
    path('api/blast-tasks/<uuid:task_id>/download/identity-plot/', blast_download_identity_plot, name='blast_download_identity_plot'),
    path('api/blast-tasks/<uuid:task_id>/download/all-visualizations/', blast_download_all_visualizations, name='blast_download_all_visualizations'),
    # --- 新增路由：结束 ---

    # 新增：基因表达分析功能相关URL
    path('gene-expression-analysis/', gene_expression_analysis.gene_analysis_page, name='gene_expression_analysis_page'),
    path('api/gene-tasks/', gene_expression_analysis.create_gene_task, name='create_gene_task'),
    path('api/gene-tasks/<uuid:task_id>/upload/', gene_expression_analysis.upload_gene_file, name='upload_gene_file'),
    path('api/gene-tasks/<uuid:task_id>/run/', gene_expression_analysis.run_gene_analysis, name='run_gene_analysis'),
    path('api/gene-tasks/<uuid:task_id>/status/', gene_expression_analysis.get_gene_task_status, name='get_gene_analysis_status'),
    path('api/gene-tasks/<uuid:task_id>/results/', gene_expression_analysis.get_gene_task_results, name='get_gene_analysis_results'),
    path('api/gene-results/<uuid:result_id>/download/', gene_expression_analysis.download_gene_result, name='download_gene_result'),
    path('api/gene-tasks/<uuid:task_id>/download-all/', gene_expression_analysis.download_all_gene_results, name='download_all_gene_results'),
    
    # 新增：DNA甲基化分析功能相关URL
    path('dna-methylation-analysis/', dna_methylation_analysis_page, name='dna_methylation_analysis_page'),
    path('api/dna-tasks/', dna_create_task, name='dna_create_task'),
    path('api/dna-tasks/<uuid:task_id>/upload/', dna_upload_file, name='dna_upload_file'),
    path('api/dna-tasks/<uuid:task_id>/run/', dna_run_analysis, name='dna_run_analysis'),
    path('api/dna-tasks/<uuid:task_id>/status/', dna_get_task_status, name='dna_get_task_status'),
    path('api/dna-tasks/<uuid:task_id>/results/', dna_get_task_results, name='dna_get_task_results'),
    path('api/dna-results/<uuid:result_id>/download/', dna_download_result, name='dna_download_result'),
    path('api/dna-tasks/<uuid:task_id>/download-all/', dna_download_all_results, name='dna_download_all_results'),
    path('api/dna-tasks/<uuid:task_id>/plots/<str:plot_type>/', dna_get_visualization_plot, name='dna_get_visualization_plot'),
    path('api/dna-tasks/<uuid:task_id>/download-by-type/', dna_download_result_by_type, name='dna_download_result_by_type'),
    path('api/dna-test/<uuid:task_id>/', test_dna_url, name='dna_test_url'),
    
    # ===== 转录组分析相关的URL =====
    # 转录组分析页面
    path('transcriptome/', transcriptome_analysis_page, name='transcriptome_page'),
    
    # API接口 - 统一使用api前缀和transcriptome-tasks格式
    # 任务管理
    path('api/transcriptome-tasks/', create_transcriptome_task, name='create_transcriptome_task'),
    path('api/transcriptome-tasks/<uuid:task_id>/cancel/', cancel_transcriptome_task, name='cancel_transcriptome_task'),
    path('api/transcriptome-tasks/list/', list_user_transcriptome_tasks, name='list_user_transcriptome_tasks'),
    
    # 文件上传
    path('api/transcriptome-tasks/<uuid:task_id>/upload/', upload_transcriptome_file, name='upload_transcriptome_file'),
    
    # 分析执行
    path('api/transcriptome-tasks/<uuid:task_id>/run/', run_transcriptome_analysis, name='run_transcriptome_analysis'),
    
    # 状态查询
    path('api/transcriptome-tasks/<uuid:task_id>/status/', get_transcriptome_task_status, name='get_transcriptome_task_status'),
    
    # 结果获取
    path('api/transcriptome-tasks/<uuid:task_id>/results/', get_transcriptome_task_results, name='get_transcriptome_task_results'),
    
    # 可视化图表
    path('api/transcriptome-tasks/<uuid:task_id>/plots/<str:plot_type>/', get_transcriptome_visualization_plot, name='get_transcriptome_plot'),
    
    # 结果下载
    path('api/transcriptome-tasks/<uuid:task_id>/download/', download_transcriptome_result_by_type, name='download_transcriptome_result_by_type'),
    path('api/transcriptome-tasks/<uuid:task_id>/download-all/', download_all_transcriptome_results, name='download_all_transcriptome_results'),
    # ===== 转录组分析URL结束 =====
    
    # KEGG分析页面
    path('kegg-pathway-analysis/', kegg_analysis_page, name='kegg_analysis_page'),
    # API接口 - KEGG分析
    # 任务管理
    path('api/kegg-tasks/', create_kegg_task, name='create_kegg_task'),
    path('api/kegg-tasks/<uuid:task_id>/cancel/', cancel_kegg_task, name='cancel_kegg_task'),
    path('api/kegg-tasks/list/', list_user_kegg_tasks, name='list_user_kegg_tasks'),
    # 文件上传和验证
    path('api/kegg/validate/', validate_kegg_gene_list, name='validate_kegg_gene_list'),
    path('api/kegg/previous-results/', get_kegg_previous_results, name='get_kegg_previous_results'),
    path('api/kegg/save-params/', save_kegg_parameters, name='save_kegg_parameters'),
    # 文件上传和任务执行
    path('api/kegg-tasks/<uuid:task_id>/upload/', upload_kegg_file, name='upload_kegg_file'),
    path('api/kegg-tasks/<uuid:task_id>/run/', run_kegg_analysis, name='run_kegg_analysis'),
    # 状态查询和结果获取
    path('api/kegg-tasks/<uuid:task_id>/status/', get_kegg_task_status, name='get_kegg_task_status'),
    path('api/kegg-tasks/<uuid:task_id>/results/', get_kegg_task_results, name='get_kegg_task_results'),
    path('api/kegg/pathway/<str:pathway_id>/', get_kegg_pathway_detail, name='get_kegg_pathway_detail'),
    # 可视化图表
    path('api/kegg-tasks/<uuid:task_id>/plots/<str:plot_type>/', get_kegg_visualization_plot, name='get_kegg_plot'),
    # 结果下载
    path('api/kegg-tasks/<uuid:task_id>/download/', download_kegg_result_by_type, name='download_kegg_result'),
    path('api/kegg-tasks/<uuid:task_id>/download-all/', download_all_kegg_results, name='download_all_kegg_results'),
    path('api/kegg-tasks/<uuid:task_id>/download-report/', download_kegg_analysis_report, name='download_kegg_report'),
    path('api/kegg/preview-mapping/<uuid:task_id>/', preview_kegg_mapping, name='preview_kegg_mapping'),
    path('api/kegg/database-stats/', get_kegg_database_stats, name='kegg_database_stats'),


    # ai功能
    path('api/deepseek/chat', csrf_exempt(views.UnifiedDeepSeekView.as_view()), name='deepseek_chat'),
    path('generate-report/', views.generate_hormone_report, name='generate_hormone_report'),
    path('api/extract_entity/', views.extract_entity_for_pdf_view, name='extract_entity'),
    # 基因预测
    path('gene-prediction/', views.gene_prediction, name='gene_prediction'),
    path('download-pdf/', views.download_pdf, name='download_pdf'),
    
    
    
     # GNN预测相关URL
    path('gnn/', gnn_views.gnn_dashboard, name='gnn_dashboard'),
    path('api/gnn/predictions/', gnn_views.gnn_get_predictions, name='gnn_get_predictions'),
    path('api/gnn/predictions/<uuid:prediction_id>/', gnn_views.gnn_get_prediction_detail, name='gnn_prediction_detail'),
    path('api/gnn/network-graph/', gnn_views.gnn_get_network_graph, name='gnn_network_graph'),
    path('api/gnn/download/csv/', gnn_views.gnn_download_predictions_csv, name='gnn_download_csv'),
    path('api/gnn/download/viz/<str:viz_type>/', gnn_views.gnn_download_visualization, name='gnn_download_viz'),
    
    # 跨物种比较相关URL
    path('phylogeny/', phylogeny_views.phylogeny_dashboard, name='phylogeny_dashboard'),
    path('api/phylogeny/trees/', phylogeny_views.phylogeny_list_trees, name='phylogeny_list_trees'),
    path('api/phylogeny/trees/<int:tree_id>/', phylogeny_views.phylogeny_get_tree, name='phylogeny_get_tree'),
    path('api/phylogeny/trees/<int:tree_id>/download/', phylogeny_views.phylogeny_download_tree, name='phylogeny_download_tree'),
    path('api/phylogeny/build-tree/', phylogeny_views.phylogeny_build_tree, name='phylogeny_build_tree'),
    path('api/phylogeny/conservation/', phylogeny_views.phylogeny_get_conservation, name='phylogeny_conservation'),
    path('api/phylogeny/families/', phylogeny_views.phylogeny_get_families, name='phylogeny_families'),
    path('api/phylogeny/download/conservation-csv/', phylogeny_views.phylogeny_download_conservation_csv, name='phylogeny_download_conservation_csv'),
    path('api/phylogeny/download/heatmap/', phylogeny_views.phylogeny_download_heatmap, name='phylogeny_download_heatmap'),

# See页面
    
    # 论文问答API
    path('page11/', see_view.page11, name='page11'),
    
    path('api/hormone-data/', views.api_hormone_data, name='api_hormone_data'),
    path('api/hormone-stats/', views.api_hormone_stats, name='api_hormone_stats'),
    path('api/network-data/', views.network_data, name='network_data'),
    path('api/paper-qa/', see_view.paper_qa, name='paper_qa'),
    path('api/chart-data/', views.chart_data, name='chart_data'),
    path("api/v263-sequences/", views.api_v263_sequence_index, name="api_v263_sequence_index"),
    path("api/v263-sequences/stats/", views.api_v263_sequence_stats, name="api_v263_sequence_stats"),
]
