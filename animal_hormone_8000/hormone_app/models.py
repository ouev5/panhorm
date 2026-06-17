from django.db import models
from django.contrib.auth.models import User
import uuid
from django.utils import timezone
import os
from django.conf import settings

# 确保媒体文件目录存在（在实际部署中可能更适合放在信号或启动脚本中）
def ensure_media_directories():
    os.makedirs(os.path.join(settings.MEDIA_ROOT, 'gene_uploads'), exist_ok=True)
    os.makedirs(os.path.join(settings.MEDIA_ROOT, 'gene_results'), exist_ok=True)

# 调用函数确保目录存在
ensure_media_directories()


class Hormone(models.Model):
    CATEGORY_CHOICES = [
        ('PEP', '肽类激素'),
        ('STE', '类固醇激素'),
        ('AMI', '胺类激素')
    ]
    id = models.AutoField(primary_key=True)
    receptor_type = models.CharField(max_length=255, verbose_name="受体类型", blank=True, null=True)
    receptor_name = models.CharField(max_length=255, verbose_name="受体名称")
    gene_family = models.CharField(max_length=255, verbose_name="基因家族", blank=True, null=True)
    gene_name_symbol = models.CharField(max_length=255, verbose_name="基因名称/符号", blank=True, null=True)
    ligand = models.CharField(max_length=255, verbose_name="配体", blank=True, null=True)
    ligand_type_comments = models.TextField(verbose_name="配体类型/注释", blank=True, null=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        db_table = "hormone_receptors_table"
        verbose_name = "激素受体数据"
        verbose_name_plural = verbose_name

    def __str__(self):
        return self.receptor_name


class Literature(models.Model):
    title = models.CharField(max_length=255, verbose_name='标题')
    authors = models.CharField(max_length=255, verbose_name='作者')
    year = models.IntegerField(verbose_name='发表年份')
    abstract = models.TextField(blank=True, null=True, verbose_name='摘要')
    file = models.FileField(upload_to='literature_files/', verbose_name='文献文件')
    upload_date = models.DateTimeField(default=timezone.now, verbose_name='上传时间')
    lit_type = models.CharField(max_length=50, choices=[
        ('ARTICLE', '期刊论文'),
        ('THESIS', '学位论文'),
        ('REPORT', '研究报告'),
        ('OTHER', '其他')
    ], default='ARTICLE', verbose_name='文献类型')

    class Meta:
        verbose_name = '文献'
        verbose_name_plural = '文献'
        ordering = ['-upload_date']

    def __str__(self):
        return f"{self.title} ({self.year})"


class AnalysisHistory(models.Model):
    """用户基因分析历史记录"""
    user = models.ForeignKey(User, on_delete=models.CASCADE, verbose_name="用户")
    sequence = models.TextField(verbose_name="分析序列", max_length=10000)
    results = models.JSONField(verbose_name="分析结果")
    created_at = models.DateTimeField(auto_now_add=True, verbose_name="分析时间")

    class Meta:
        verbose_name = "分析历史"
        verbose_name_plural = "分析历史"
        ordering = ['-created_at']

    def __str__(self):
        return f"分析 #{self.id} by {self.user.username}"


class AnalysisTask(models.Model):
    """异步分析任务模型"""
    TASK_TYPES = (
        ('enrichment', '富集分析'),
        ('neural', '神经网络分析'),
        ('coexpression', '激素-基因共表达网络分析'),
        ('evolution', '进化保守性分析'),
        ('timeseries', '动态时序建模'),
        ('phenotype', '表型关联富集分析'),
    )
    
    STATUS_CHOICES = (
        ('pending', '等待中'),
        ('running', '进行中'),
        ('completed', '已完成'),
        ('failed', '失败'),
    )
    
    task_id = models.CharField(max_length=36, default=uuid.uuid4, editable=False, unique=True)
    analysis_type = models.CharField(max_length=20, choices=TASK_TYPES)
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='pending')
    input_file = models.FileField(upload_to='uploads/')
    progress = models.IntegerField(default=0)
    result = models.JSONField(null=True, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    
    def __str__(self):
        return f"{self.analysis_type} - {self.task_id}"    


class UserProfile(models.Model):
    """用户扩展资料"""
    user = models.OneToOneField(User, on_delete=models.CASCADE)
    institution = models.CharField(max_length=255, blank=True, null=True, verbose_name="机构")
    position = models.CharField(max_length=100, blank=True, null=True, verbose_name="职位")
    avatar = models.ImageField(upload_to='avatars/', blank=True, null=True, verbose_name="头像")
    bio = models.TextField(blank=True, null=True, verbose_name="个人简介")

    class Meta:
        verbose_name = "用户资料"
        verbose_name_plural = "用户资料"

    def __str__(self):
        return self.user.username    


class HormoneRawData(models.Model):
    id = models.IntegerField(primary_key=True, verbose_name="ID")
    hormone = models.TextField(verbose_name="激素名称/结构")
    receptors_in_article = models.CharField(max_length=10000, verbose_name="文章中的受体", blank=True, null=True)
    gene_name = models.CharField(max_length=255, verbose_name="基因名称", blank=True, null=True)
    gene_sequence = models.TextField(verbose_name="基因序列", blank=True, null=True)
    pmid = models.CharField(max_length=255, verbose_name="PMID", blank=True, null=True)
    created_at = models.DateTimeField(auto_now_add=True, verbose_name="创建时间")
    updated_at = models.DateTimeField(auto_now=True, verbose_name="更新时间")

    class Meta:
        db_table = "hormone_raw_data"  # 数据库表名，可自定义
        verbose_name = "激素原始数据"
        verbose_name_plural = "激素原始数据"

    def __str__(self):
        return f"{self.id} - {self.hormone}"


class HormoneRelatedGene(models.Model):
    id = models.IntegerField(primary_key=True)  # 主键ID，对应数据库NOT NULL
    hormone_name = models.CharField(max_length=255, null=True)  # 激素名称，允许为NULL
    related_genes = models.CharField(max_length=255, null=True)  # 相关基因，允许为NULL
    pmid = models.BigIntegerField(null=True)  # PubMed ID，对应数据库bigint类型，允许为NULL
    gene_sequence = models.TextField(null=True)  # 基因序列，允许为NULL
    related_diseases = models.TextField(null=True)  # 相关疾病，对应数据库text类型，允许为NULL
    # 新增的三个字段，根据你的数据特征选择合适的字段类型
    DO_Match_disease_names = models.CharField(max_length=255, null=True, verbose_name="DO匹配疾病名称")
    DO_ID = models.CharField(max_length=50, null=True, verbose_name="DO编号")  # DOID格式固定，长度无需过长
    DO_Standardized_Terminology = models.TextField(null=True, verbose_name="DO标准化术语")  # 术语可能较长，用TextField
    hormone_accession = models.CharField(max_length=50, null=True, verbose_name="激素UniProt ID")
    organism = models.CharField(max_length=255, null=True, verbose_name="物种")
    regulation_type = models.CharField(max_length=50, null=True, verbose_name="调控类型")
    pmid_gene_disease = models.CharField(max_length=50, null=True, verbose_name="基因-疾病PMID")
    relationship_gene_disease = models.CharField(max_length=100, null=True, verbose_name="基因-疾病关系类型")
    has_pmid_gene_disease = models.CharField(max_length=10, null=True, verbose_name="是否有基因-疾病PMID")
    doid_standardized_name = models.CharField(max_length=500, null=True, verbose_name="DO标准疾病名称")
    evidence_source = models.CharField(max_length=100, null=True, verbose_name="证据来源")
    # 新增基因序列字段（2026-04-28）
    gene_description = models.TextField(null=True, verbose_name="基因描述")
    chromosome = models.CharField(max_length=50, null=True, verbose_name="染色体")
    map_location = models.CharField(max_length=100, null=True, verbose_name="图谱位置")
    gene_synonyms = models.TextField(null=True, verbose_name="基因别名")
    chr_accession = models.CharField(max_length=100, null=True, verbose_name="染色体编号")
    protein_accession = models.CharField(max_length=100, null=True, verbose_name="蛋白质编号")
    protein_sequence = models.TextField(null=True, verbose_name="蛋白质序列")
    mrna_accession = models.CharField(max_length=100, null=True, verbose_name="mRNA编号")
    mrna_sequence = models.TextField(null=True, verbose_name="mRNA序列")

    class Meta:
        db_table = 'hormone_data'  # 数据库表名
        verbose_name = '激素相关基因'
        verbose_name_plural = '激素相关基因'

    def __str__(self):
        return f"{self.hormone_name} - {self.related_genes}"


class HormoneReceptorInfo(models.Model):
    # 主键字段移除null=True，保留db_column和primary_key=True
    id = models.TextField(db_column='ID', primary_key=True, blank=True)  # 移除null=True
    pubchem_id = models.TextField(db_column='Pubchem ID', blank=True, null=True)
    hormone_name = models.TextField(db_column='hormone name', blank=True, null=True)
    receptor_uniprot_id = models.TextField(db_column='receptor uniProt ID', blank=True, null=True)
    receptor_coding_genes = models.TextField(db_column='receptor coding genes', blank=True, null=True)
    receptor_name = models.TextField(db_column='receptor name', blank=True, null=True)
    receptor_species_name = models.TextField(db_column='receptor species name', blank=True, null=True)
    receptor_coding_genes_sequence = models.TextField(db_column='receptor coding genes sequence', blank=True, null=True)

    class Meta:
        db_table = 'hormone_feitailei'
        verbose_name = '非肽类激素受体信息'
        verbose_name_plural = '非肽类激素受体信息'

    def __str__(self):
        return f"{self.hormone_name or '未知激素'} - {self.receptor_name or '未知受体'}"

class HormoneReceptorFull(models.Model):
    # 1. 字段名通过db_column映射到数据库的带空格列名
    # 2. 数据库类型是text，模型用TextField匹配
    # 3. 数据库允许空（Null=YES），模型添加blank=True, null=True
    id = models.CharField(max_length=255, db_column='ID', primary_key=True)
    hormone_uniprot_id = models.TextField(db_column='hormone uniProt ID', blank=True, null=True)
    hormone_name = models.TextField(db_column='hormone name', blank=True, null=True)  # 映射“hormone name”列
    hormone_species_name = models.TextField(db_column='hormone species name', blank=True, null=True)
    hormone_coding_genes = models.TextField(db_column='hormone coding genes', blank=True, null=True)
    hormone_coding_genes_sequence = models.TextField(db_column='hormone coding genes sequence', blank=True, null=True)
    receptor_uniprot_id = models.TextField(db_column='receptor uniProt ID', blank=True, null=True)
    receptor_name = models.TextField(db_column='receptor name', blank=True, null=True)
    receptor_species_name = models.TextField(db_column='receptor species name', blank=True, null=True)
    receptor_coding_genes = models.TextField(db_column='receptor coding genes', blank=True, null=True)
    receptor_coding_genes_sequence = models.TextField(db_column='receptor coding genes sequence', blank=True, null=True)

    class Meta:
        db_table = 'hormone_tailei'  # 与数据库表名一致
        verbose_name = '肽类激素受体信息'
        verbose_name_plural = '肽类激素受体信息'

    def __str__(self):
        return f"{self.hormone_name or '未知激素'} - {self.receptor_name or '未知受体'}"


# 基因分析相关模型
class GeneAnalysisTask(models.Model):
    """基因表达差异分析任务模型"""
    TASK_STATUS = (
        ('pending', '等待中'),
        ('running', '运行中'),
        ('completed', '已完成'),
        ('failed', '失败')
    )
    
    task_id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    user = models.ForeignKey(User, on_delete=models.CASCADE, null=True, blank=True)
    task_name = models.CharField(max_length=200, default="未命名基因分析")
    status = models.CharField(max_length=20, choices=TASK_STATUS, default='pending')
    progress = models.IntegerField(default=0)
    created_at = models.DateTimeField(auto_now_add=True)
    completed_at = models.DateTimeField(null=True, blank=True)
    error_message = models.TextField(null=True, blank=True)
    
    # 分析参数
    data_params = models.JSONField(default=dict, blank=True)
    normalization_params = models.JSONField(default=dict, blank=True)
    de_params = models.JSONField(default=dict, blank=True)
    visualization_params = models.JSONField(default=dict, blank=True)
    pathway_params = models.JSONField(default=dict, blank=True)
    
    def __str__(self):
        return f"{self.task_name} ({self.task_id})"


class GeneUploadedFile(models.Model):
    """基因分析上传文件模型"""
    FILE_TYPES = (
        ('expression-data', '表达数据'),
        ('sample-info', '样本信息'),
        ('gene-annotations', '基因注释'),
        ('pathway-data', '通路数据'),
        ('contrast-definitions', '对比定义')
    )
    
    file_id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    task = models.ForeignKey(GeneAnalysisTask, related_name='uploaded_files', on_delete=models.CASCADE)
    file_type = models.CharField(max_length=20, choices=FILE_TYPES)
    file = models.FileField(upload_to='gene_uploads/')
    filename = models.CharField(max_length=255)
    file_size = models.IntegerField(help_text="文件大小(字节)")
    uploaded_at = models.DateTimeField(auto_now_add=True)
    
    def __str__(self):
        return f"{self.get_file_type_display()}: {self.filename}"


class GeneAnalysisResult(models.Model):
    """基因分析结果模型"""
    RESULT_TYPES = (
        ('normalized-data', '标准化数据'),
        ('de-results', '差异表达结果'),
        ('volcano-plot', '火山图'),
        ('heatmap', '热图'),
        ('pathway-analysis', '通路分析结果'),
        ('report', '分析报告')
    )
    
    result_id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    task = models.ForeignKey(GeneAnalysisTask, related_name='results', on_delete=models.CASCADE)
    result_type = models.CharField(max_length=30, choices=RESULT_TYPES)
    file = models.FileField(upload_to='gene_results/')
    filename = models.CharField(max_length=255)
    created_at = models.DateTimeField(auto_now_add=True)
    
    def __str__(self):
        return f"{self.get_result_type_display()}: {self.filename}"

class AnalysisJob(models.Model):
    """分析任务主模型"""
    JOB_STATUS_CHOICES = [
        ('pending', 'Pending'),
        ('running', 'Running'),
        ('completed', 'Completed'),
        ('failed', 'Failed'),
    ]
    
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    user = models.ForeignKey(User, on_delete=models.CASCADE, related_name='analysis_jobs')
    job_name = models.CharField(max_length=255, default="Untitled Analysis")
    status = models.CharField(max_length=20, choices=JOB_STATUS_CHOICES, default='pending')
    created_at = models.DateTimeField(default=timezone.now)
    started_at = models.DateTimeField(null=True, blank=True)
    completed_at = models.DateTimeField(null=True, blank=True)
    progress = models.IntegerField(default=0)  # 0-100
    log = models.TextField(blank=True)
    
    def __str__(self):
        return f"{self.job_name} ({self.user.username})"


class UploadedFile(models.Model):
    """上传的文件模型"""
    FILE_TYPE_CHOICES = [
        ('methyl_data', 'Methylation Data'),
        ('metadata', 'Sample Metadata'),
        ('regions', 'Genomic Regions'),
        ('expression', 'Expression Data'),
        ('custom_genesets', 'Custom Gene Sets'),
    ]
    
    job = models.ForeignKey(AnalysisJob, on_delete=models.CASCADE, related_name='uploaded_files')
    file_type = models.CharField(max_length=20, choices=FILE_TYPE_CHOICES)
    file = models.FileField(upload_to='uploads/%Y/%m/%d/')
    filename = models.CharField(max_length=255)
    uploaded_at = models.DateTimeField(default=timezone.now)
    file_size = models.IntegerField(help_text="File size in bytes")
    
    def __str__(self):
        return f"{self.get_file_type_display()} - {self.filename}"


class AnalysisParameter(models.Model):
    """分析参数模型"""
    PARAMETER_GROUP_CHOICES = [
        ('quality_control', 'Quality Control'),
        ('preprocessing', 'Preprocessing'),
        ('differential', 'Differential Methylation'),
        ('region_analysis', 'Region Analysis'),
        ('correlation', 'Correlation'),
        ('enrichment', 'Enrichment'),
    ]
    
    job = models.ForeignKey(AnalysisJob, on_delete=models.CASCADE, related_name='parameters')
    group = models.CharField(max_length=20, choices=PARAMETER_GROUP_CHOICES)
    name = models.CharField(max_length=100)
    value = models.TextField()
    
    class Meta:
        unique_together = ('job', 'group', 'name')
    
    def __str__(self):
        return f"{self.job.job_name} - {self.get_group_display()} - {self.name}"


class AnalysisResult(models.Model):
    """分析结果模型"""
    RESULT_TYPE_CHOICES = [
        ('dmp', 'Differentially Methylated Positions'),
        ('dmr', 'Differentially Methylated Regions'),
        ('enrichment', 'Enrichment Results'),
        ('methyl_values', 'Methylation Values'),
        ('visualization', 'Visualization'),
        ('summary', 'Summary Statistics'),
    ]
    
    job = models.ForeignKey(AnalysisJob, on_delete=models.CASCADE, related_name='results')
    result_type = models.CharField(max_length=20, choices=RESULT_TYPE_CHOICES)
    file = models.FileField(upload_to='results/%Y/%m/%d/')
    filename = models.CharField(max_length=255)
    created_at = models.DateTimeField(default=timezone.now)
    description = models.TextField(blank=True)
    
    def __str__(self):
        return f"{self.get_result_type_display()} - {self.filename}"

class BlastTask(models.Model):
    STATUS_CHOICES = [
        ('created', 'Created'),
        ('pending', 'Pending'),
        ('running', 'Running'),
        ('completed', 'Completed'),
        ('failed', 'Failed'),
    ]

    id = models.CharField(max_length=255, primary_key=True)  # UUID
    name = models.CharField(max_length=255)
    description = models.TextField(blank=True)
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='created')
    sequence_type = models.CharField(max_length=20)  # 'protein' or 'nucleotide'
    error_message = models.TextField(blank=True, null=True, verbose_name='错误信息')
    
    # Parameters
    sequence_file = models.CharField(max_length=500, blank=True, null=True)
    
    # 修正：增加字段长度
    database_requested = models.CharField(max_length=100, blank=True, null=True)  # 增加长度
    database_actual = models.CharField(max_length=255, blank=True, null=True)     # 大幅增加长度
    
    program = models.CharField(max_length=20, blank=True, null=True)
    evalue = models.CharField(max_length=20, blank=True, null=True)
    max_hits = models.CharField(max_length=20, blank=True, null=True)
    
    # Result files
    result_file = models.CharField(max_length=500, blank=True, null=True)
    csv_result_file = models.CharField(max_length=500, blank=True, null=True)
    fasta_hits_file = models.CharField(max_length=500, blank=True, null=True)
    xml_report_file = models.CharField(max_length=500, blank=True, null=True)

    created_at = models.DateTimeField(default=timezone.now)
    completed_at = models.DateTimeField(blank=True, null=True)

    def __str__(self):
        return f"BlastTask({self.name}, {self.status})"
    
    # 可以添加一些有用的方法
    def get_duration(self):
        """获取任务执行时长"""
        if self.completed_at and self.created_at:
            return self.completed_at - self.created_at
        return None
    
    def is_completed(self):
        return self.status == 'completed'
    
    def is_failed(self):
        return self.status == 'failed'
    
    def is_running(self):
        return self.status == 'running'

class DNAMethylationTask(models.Model):
    task_id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    status = models.CharField(max_length=50, default='created')
    progress = models.IntegerField(default=0)
    created_at = models.DateTimeField(auto_now_add=True)
    completed_at = models.DateTimeField(null=True, blank=True)
    params = models.JSONField(default=dict, blank=True)
    files = models.JSONField(default=list, blank=True)
    results = models.JSONField(null=True, blank=True)
    error = models.TextField(blank=True)
    
    class Meta:
        db_table = 'dna_methylation_tasks'
        
class TranscriptomeAnalysisTask(models.Model):
    """转录组分析任务模型"""
    
    TASK_STATUS = (
        ('pending', '待处理'),
        ('running', '运行中'),
        ('completed', '已完成'),
        ('failed', '失败'),
        ('cancelled', '已取消')
    )
    
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    user = models.ForeignKey(User, on_delete=models.CASCADE, null=True, blank=True)
    username = models.CharField(max_length=100, default='anonymous')
    
    status = models.CharField(max_length=20, choices=TASK_STATUS, default='pending')
    parameters = models.JSONField(default=dict, blank=True)  # 分析参数
    files = models.JSONField(default=dict, blank=True)  # 文件信息
    analysis_type = models.CharField(max_length=50, default='transcriptome')
    
    progress = models.IntegerField(default=0)  # 进度百分比
    current_step = models.CharField(max_length=100, default='Initializing')
    message = models.TextField(blank=True, null=True)
    error = models.TextField(blank=True, null=True)
    results = models.JSONField(default=dict, blank=True, null=True)  # 分析结果
    
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    started_at = models.DateTimeField(null=True, blank=True)
    completed_at = models.DateTimeField(null=True, blank=True)
    
    class Meta:
        db_table = 'transcriptome_analysis_tasks'
        ordering = ['-created_at']
        indexes = [
            models.Index(fields=['user', 'status']),
            models.Index(fields=['created_at']),
        ]
    
    def __str__(self):
        return f"Transcriptome Task {self.id} - {self.status}"


class TranscriptomeAnalysisFile(models.Model):
    """转录组分析文件模型"""
    
    FILE_TYPES = (
        ('fastq-files', 'FASTQ文件'),
        ('sample-metadata', '样本元数据'),
        ('reference-genome', '参考基因组'),
        ('gene-list', '基因列表'),
        ('custom-genesets', '自定义基因集')
    )
    
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    task = models.ForeignKey(TranscriptomeAnalysisTask, on_delete=models.CASCADE, related_name='files_relation')
    
    file_name = models.CharField(max_length=255)
    file_path = models.TextField()  # 服务器上的文件路径
    file_type = models.CharField(max_length=50, choices=FILE_TYPES)
    file_size = models.BigIntegerField()  # 文件大小（字节）
    
    uploaded_at = models.DateTimeField(auto_now_add=True)
    
    class Meta:
        db_table = 'transcriptome_analysis_files'
        ordering = ['uploaded_at']
        indexes = [
            models.Index(fields=['task', 'file_type']),
        ]
    
    def __str__(self):
        return f"{self.file_name} ({self.file_type})"


class TranscriptomeAnalysisResult(models.Model):
    """转录组分析结果模型"""
    
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    task = models.OneToOneField(TranscriptomeAnalysisTask, on_delete=models.CASCADE, related_name='result_relation')
    
    summary = models.JSONField(default=dict, blank=True)  # 分析摘要
    differential_expression = models.JSONField(default=dict, blank=True)  # 差异表达结果
    enrichment = models.JSONField(default=dict, blank=True)  # 富集分析结果
    file_paths = models.JSONField(default=dict, blank=True)  # 结果文件路径
    
    created_at = models.DateTimeField(auto_now_add=True)
    
    class Meta:
        db_table = 'transcriptome_analysis_results'
        ordering = ['-created_at']
    
    def __str__(self):
        return f"Results for Task {self.task.id}"
        
class KeggAnalysisTask(models.Model):
    STATUS_CHOICES = [
        ('created', '已创建'),
        ('uploaded', '已上传文件'),
        ('running', '分析中'),
        ('completed', '已完成'),
        ('failed', '失败'),
        ('cancelled', '已取消'),
    ]
    
    task_id = models.UUIDField(default=uuid.uuid4, editable=False, unique=True)
    client_id = models.CharField(max_length=255)  # 客户端标识
    client_ip = models.GenericIPAddressField()
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='created')
    parameters = models.JSONField(default=dict, blank=True)
    genes = models.JSONField(default=list, blank=True)  # 存储基因列表
    gene_count = models.IntegerField(default=0)
    file_name = models.CharField(max_length=255, blank=True, null=True)
    result_file = models.FileField(upload_to='kegg_results/', blank=True, null=True)
    results = models.JSONField(default=dict, blank=True)  # 存储分析结果
    progress = models.IntegerField(default=0)
    message = models.TextField(blank=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    
    class Meta:
        ordering = ['-created_at']
    
    def __str__(self):
        return f"KEGG Task {self.task_id} ({self.status})"
        
from django.db import models
from django.contrib.postgres.fields import ArrayField  # 如果用PostgreSQL
import uuid
from django.core.exceptions import ValidationError

class KeggOrganism(models.Model):
    """KEGG生物体"""
    code = models.CharField(max_length=10, primary_key=True, verbose_name="KEGG代码")
    name = models.CharField(max_length=200, verbose_name="生物体名称")
    scientific_name = models.CharField(max_length=200, blank=True, null=True, verbose_name="学名")  # 增加null=True
    description = models.TextField(blank=True, null=True, verbose_name="描述")  # 增加null=True
    taxonomy = models.TextField(blank=True, null=True, verbose_name="分类信息")  # 增加null=True
    
    class Meta:
        verbose_name = "KEGG生物体"
        verbose_name_plural = "KEGG生物体"
        indexes = [models.Index(fields=['code'])]  # 新增索引优化查询
    
    def __str__(self):
        return f"{self.code} - {self.name}"

class KeggPathway(models.Model):
    """KEGG通路"""
    pathway_id = models.CharField(max_length=20, primary_key=True, verbose_name="通路ID")
    name = models.CharField(max_length=500, verbose_name="通路名称")
    description = models.TextField(blank=True, null=True, verbose_name="通路描述")  # 增加null=True
    
    # 分类信息
    CATEGORY_CHOICES = [
        ('metabolism', '代谢'),
        ('genetic', '遗传信息处理'),
        ('environmental', '环境信息处理'),
        ('cellular', '细胞过程'),
        ('organismal', '生物体系统'),
        ('human', '人类疾病'),
        ('drug', '药物开发'),
    ]
    category = models.CharField(max_length=50, choices=CATEGORY_CHOICES, blank=True, null=True, verbose_name="分类")  # 增加空值处理
    subcategory = models.CharField(max_length=100, blank=True, null=True, verbose_name="子分类")  # 增加null=True
    
    # 统计信息
    gene_count = models.IntegerField(default=0, verbose_name="基因数量")
    compound_count = models.IntegerField(default=0, verbose_name="化合物数量")
    reaction_count = models.IntegerField(default=0, verbose_name="反应数量")
    
    # 链接信息
    kegg_url = models.URLField(blank=True, null=True, verbose_name="KEGG链接")  # 增加null=True
    image_url = models.URLField(blank=True, null=True, verbose_name="通路图链接")  # 增加null=True
    
    # 时间戳
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    
    class Meta:
        verbose_name = "KEGG通路"
        verbose_name_plural = "KEGG通路"
        indexes = [
            models.Index(fields=['category']),
            models.Index(fields=['pathway_id']),
            models.Index(fields=['name']),  # 新增名称索引，优化模糊查询
        ]
    
    def __str__(self):
        return f"{self.pathway_id} - {self.name}"

class KeggGene(models.Model):
    """KEGG基因（核心优化：解决symbol字段长度超限问题）"""
    gene_id = models.CharField(max_length=50, primary_key=True, verbose_name="KEGG基因ID")
    
    # 关键修改：扩大symbol字段长度，增加空值处理，适配KEGG超长基因符号
    symbol = models.CharField(
        max_length=200,  # 从50→200，解决Data too long错误
        db_index=True, 
        blank=True, 
        null=True,  # 允许空值，避免数据缺失时插入失败
        verbose_name="基因符号"
    )
    
    # 优化：基因名称改用TextField，适配超长名称（KEGG部分基因名称超过500字符）
    name = models.TextField(blank=True, null=True, verbose_name="基因名称")
    
    # 外部数据库ID（增加null=True，适配数据缺失场景）
    entrez_id = models.CharField(max_length=20, blank=True, null=True, db_index=True, verbose_name="Entrez ID")
    uniprot_id = models.CharField(max_length=20, blank=True, null=True, db_index=True, verbose_name="UniProt ID")
    ensembl_id = models.CharField(max_length=50, blank=True, null=True, db_index=True, verbose_name="Ensembl ID")
    
    # 所属生物体（修改删除策略为PROTECT，防止误删生物体导致基因数据丢失）
    organism = models.ForeignKey(
        'KeggOrganism',  # 使用字符串引用避免循环导入问题
        on_delete=models.PROTECT,  # 替换CASCADE为PROTECT，更安全
        related_name='genes', 
        verbose_name="生物体"
    )
    
    # 基因信息（增加null=True，适配缺失的位置数据）
    chromosome = models.CharField(max_length=50, blank=True, null=True, verbose_name="染色体")  # 长度从20→50，适配特殊命名
    start_position = models.BigIntegerField(null=True, blank=True, verbose_name="起始位置")
    end_position = models.BigIntegerField(null=True, blank=True, verbose_name="终止位置")
    strand = models.CharField(max_length=1, blank=True, null=True, verbose_name="链方向")  # 增加null=True
    
    # 功能描述
    function = models.TextField(blank=True, null=True, verbose_name="功能描述")  # 增加null=True
    pathway_count = models.IntegerField(default=0, verbose_name="通路数量")
    
    # 时间戳
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    
    class Meta:
        verbose_name = "KEGG基因"
        verbose_name_plural = "KEGG基因"
        indexes = [
            models.Index(fields=['symbol']),
            models.Index(fields=['organism']),
            models.Index(fields=['entrez_id']),
            models.Index(fields=['organism', 'symbol']),  # 新增复合索引，优化常用查询
        ]
        # 移除MySQL不支持的带条件唯一约束，改为代码层验证
    
    def clean(self):
        """
        自定义数据验证：实现同一生物体下非空symbol唯一的业务规则
        会在调用full_clean()/表单验证/后台保存时自动触发
        """
        super().clean()
        
        # 仅当symbol非空且非空字符串时，检查唯一性
        if self.symbol and self.symbol.strip():
            # 构造查询条件：同一生物体 + 相同symbol
            query = KeggGene.objects.filter(
                organism=self.organism,
                symbol=self.symbol.strip()
            )
            # 如果是更新操作，排除当前记录本身
            if self.pk:
                query = query.exclude(pk=self.pk)
            
            # 检查是否存在重复
            if query.exists():
                raise ValidationError(
                    f"生物体【{self.organism}】下已存在基因符号【{self.symbol.strip()}】，无法重复创建/更新"
                )
    
    def save(self, *args, **kwargs):
        """
        重写保存方法：强制触发数据验证，确保约束生效
        覆盖所有保存场景（直接save/批量保存/后台保存）
        """
        # 强制触发clean验证，有重复则抛出ValidationError
        self.full_clean()
        
        # 可选优化：统一处理symbol空值（空字符串转为None）
        if self.symbol == "":
            self.symbol = None
        
        # 执行原生保存逻辑
        super().save(*args, **kwargs)
    
    def __str__(self):
        # 优化：处理symbol为空的情况，避免AttributeError
        symbol = self.symbol if self.symbol else "未知"
        return f"{self.gene_id} - {symbol}"

class KeggGenePathway(models.Model):
    """基因-通路关联表（优化空值和索引）"""
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    gene = models.ForeignKey(
        KeggGene, 
        on_delete=models.CASCADE, 
        related_name='pathway_relations', 
        verbose_name="基因"
    )
    pathway = models.ForeignKey(
        KeggPathway, 
        on_delete=models.CASCADE, 
        related_name='gene_relations', 
        verbose_name="通路"
    )
    
    # 关联信息（增加空值处理）
    relation_type = models.CharField(max_length=50, default='member', blank=True, null=True, verbose_name="关联类型")
    evidence = models.CharField(max_length=100, blank=True, null=True, verbose_name="证据")
    
    # 时间戳
    created_at = models.DateTimeField(auto_now_add=True)
    
    class Meta:
        verbose_name = "基因-通路关联"
        verbose_name_plural = "基因-通路关联"
        unique_together = [['gene', 'pathway']]
        indexes = [
            models.Index(fields=['gene', 'pathway']),
            models.Index(fields=['relation_type']),  # 新增关联类型索引
        ]
    
    def __str__(self):
        # 优化：处理基因symbol为空的情况
        gene_symbol = self.gene.symbol if self.gene.symbol else self.gene.gene_id
        return f"{gene_symbol} -> {self.pathway.pathway_id}"

class KeggPathwayHierarchy(models.Model):
    """通路层级关系（增加空值和索引优化）"""
    parent = models.ForeignKey(
        KeggPathway, 
        on_delete=models.CASCADE, 
        related_name='children', 
        verbose_name="父通路"
    )
    child = models.ForeignKey(
        KeggPathway, 
        on_delete=models.CASCADE, 
        related_name='parents', 
        verbose_name="子通路"
    )
    level = models.IntegerField(default=0, blank=True, null=True, verbose_name="层级")  # 增加空值处理
    
    class Meta:
        verbose_name = "通路层级"
        verbose_name_plural = "通路层级"
        unique_together = [['parent', 'child']]
        indexes = [
            models.Index(fields=['parent', 'level']),  # 新增复合索引，优化层级查询
            models.Index(fields=['child', 'level']),
        ]
    
    def __str__(self):
        return f"{self.parent.pathway_id} -> {self.child.pathway_id} (Level: {self.level})"
        
# ==================== 新增：GNN预测模块 ====================
class GNPrediction(models.Model):
    """GNN预测的激素-受体相互作用"""
    CONFIDENCE_CHOICES = [
        ('high', 'High'),
        ('medium', 'Medium'),
        ('low', 'Low'),
    ]
    
    HORMONE_TYPE_CHOICES = [
        ('peptide', '肽类激素'),
        ('non_peptide', '非肽类激素'),
    ]
    
    id = models.AutoField(primary_key=True)
    
    # 激素相关字段 - 使用 CharField 而不是直接外键
    hormone_peptide_id = models.CharField(
        max_length=255,
        verbose_name="肽类激素ID",
        null=True,
        blank=True,
        db_index=True
    )
    hormone_non_peptide_id = models.CharField(
        max_length=255,
        verbose_name="非肽类激素ID",
        null=True,
        blank=True,
        db_index=True
    )
    hormone_type = models.CharField(
        max_length=20,
        choices=HORMONE_TYPE_CHOICES,
        verbose_name="激素类型",
        null=True,
        blank=True
    )
    
    # 受体相关字段 - 使用 CharField 而不是直接外键
    receptor_peptide_id = models.CharField(
        max_length=255,
        verbose_name="肽类受体ID",
        null=True,
        blank=True,
        db_index=True
    )
    receptor_non_peptide_id = models.CharField(
        max_length=255,
        verbose_name="非肽类受体ID",
        null=True,
        blank=True,
        db_index=True
    )
    receptor_type = models.CharField(
        max_length=20,
        choices=HORMONE_TYPE_CHOICES,
        verbose_name="受体类型",
        null=True,
        blank=True
    )
    
    prediction_score = models.FloatField(help_text="Prediction score (0-1)", verbose_name="预测分数")
    confidence_level = models.CharField(max_length=10, choices=CONFIDENCE_CHOICES, verbose_name="置信度", db_index=True)
    model_version = models.CharField(max_length=50, verbose_name="模型版本", default='gnn_v1.0')
    created_at = models.DateTimeField(auto_now_add=True, verbose_name="创建时间")

    class Meta:
        db_table = 'gn_predictions'
        verbose_name = "GNN预测"
        verbose_name_plural = "GNN预测"
        indexes = [
            models.Index(fields=['prediction_score']),
            models.Index(fields=['confidence_level']),
            models.Index(fields=['hormone_peptide_id']),
            models.Index(fields=['hormone_non_peptide_id']),
            models.Index(fields=['receptor_peptide_id']),
            models.Index(fields=['receptor_non_peptide_id']),
        ]

    def get_hormone(self):
        """获取激素对象"""
        from .models import HormoneReceptorFull, HormoneReceptorInfo
        
        if self.hormone_type == 'peptide' and self.hormone_peptide_id:
            try:
                return HormoneReceptorFull.objects.get(id=self.hormone_peptide_id)
            except HormoneReceptorFull.DoesNotExist:
                return None
        elif self.hormone_type == 'non_peptide' and self.hormone_non_peptide_id:
            try:
                return HormoneReceptorInfo.objects.get(id=self.hormone_non_peptide_id)
            except HormoneReceptorInfo.DoesNotExist:
                return None
        return None
    
    def get_receptor(self):
        """获取受体对象"""
        from .models import HormoneReceptorFull, HormoneReceptorInfo
        
        if self.receptor_type == 'peptide' and self.receptor_peptide_id:
            try:
                return HormoneReceptorFull.objects.get(id=self.receptor_peptide_id)
            except HormoneReceptorFull.DoesNotExist:
                return None
        elif self.receptor_type == 'non_peptide' and self.receptor_non_peptide_id:
            try:
                return HormoneReceptorInfo.objects.get(id=self.receptor_non_peptide_id)
            except HormoneReceptorInfo.DoesNotExist:
                return None
        return None
    
    def get_hormone_name(self):
        """获取激素名称"""
        hormone = self.get_hormone()
        if hormone:
            return hormone.hormone_name if hasattr(hormone, 'hormone_name') else str(hormone.id)
        return "未知激素"
    
    def get_receptor_name(self):
        """获取受体名称"""
        receptor = self.get_receptor()
        if receptor:
            return receptor.receptor_name if hasattr(receptor, 'receptor_name') else str(receptor.id)
        return "未知受体"

    def __str__(self):
        return f"{self.get_hormone_name()} -> {self.get_receptor_name()} ({self.prediction_score:.2f})"

# ==================== 新增：跨物种比较模块 ====================
class PhylogeneticTree(models.Model):
    """激素家族系统发育树"""
    id = models.AutoField(primary_key=True)
    hormone_family = models.CharField(max_length=100, db_index=True, verbose_name="激素家族")
    tree_newick = models.TextField(verbose_name="Newick格式树")
    species_count = models.IntegerField(verbose_name="物种数量")
    created_at = models.DateTimeField(auto_now_add=True, verbose_name="创建时间")

    class Meta:
        db_table = 'phylogenetic_trees'
        ordering = ['-created_at']
        verbose_name = "系统发育树"
        verbose_name_plural = "系统发育树"

    def __str__(self):
        return f"{self.hormone_family} tree ({self.species_count} species)"


class SequenceConservation(models.Model):
    """序列保守性分析结果"""
    id = models.AutoField(primary_key=True)
    hormone = models.ForeignKey(
        'HormoneReceptorFull',  # 或 HormoneRawData，根据您的数据选择
        on_delete=models.CASCADE,
        related_name='conservation',
        verbose_name="激素"
    )
    species = models.CharField(max_length=100, verbose_name="物种")
    similarity_score = models.FloatField(verbose_name="相似度分数")
    alignment_length = models.IntegerField(verbose_name="比对长度")

    class Meta:
        db_table = 'sequence_conservation'
        unique_together = [['hormone', 'species']]
        verbose_name = "序列保守性"
        verbose_name_plural = "序列保守性"

    def __str__(self):
        hormone_name = self.hormone.hormone_name if hasattr(self.hormone, 'hormone_name') else str(self.hormone.id)
        return f"{hormone_name} - {self.species} ({self.similarity_score:.2f})"