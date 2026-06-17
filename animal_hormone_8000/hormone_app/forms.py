from django import forms
from .models import AnalysisJob, UploadedFile

class AnalysisJobForm(forms.ModelForm):
    """分析任务表单"""
    class Meta:
        model = AnalysisJob
        fields = ['job_name']
        widgets = {
            'job_name': forms.TextInput(attrs={
                'class': 'form-control',
                'placeholder': 'Enter analysis name'
            })
        }


class FileUploadForm(forms.ModelForm):
    """文件上传表单"""
    class Meta:
        model = UploadedFile
        fields = ['file_type', 'file']
        widgets = {
            'file_type': forms.HiddenInput(),
            'file': forms.FileInput(attrs={
                'class': 'form-control',
                'accept': '.bed,.bismark,.csv,.tsv,.txt,.gff,.gtf,.gmt'
            })
        }


class QualityControlForm(forms.Form):
    """质量控制参数表单"""
    minCpGsPerSample = forms.IntegerField(min_value=10000, initial=50000)
    minMeanCoverage = forms.IntegerField(min_value=1, max_value=20, initial=5)
    maxMethDeviation = forms.FloatField(min_value=0.1, max_value=0.5, step_size=0.05, initial=0.3)
    detectionPvalThresh = forms.FloatField(min_value=0.001, max_value=0.05, step_size=0.001, initial=0.01)
    minCoveragePerCpG = forms.IntegerField(min_value=1, max_value=10, initial=3)
    minSamplesWithData = forms.IntegerField(min_value=50, max_value=100, initial=70)
    excludeSnpProbes = forms.ChoiceField(choices=[
        ('yes', 'Yes, exclude all probes with SNPs'),
        ('common', 'Exclude probes with common SNPs (MAF > 0.05)'),
        ('no', 'No, keep all probes')
    ], initial='yes')
    snpDistanceThresh = forms.IntegerField(min_value=0, max_value=50, initial=5)
    batchVariable = forms.CharField(required=False, initial='batch')
    pcaComponentsBatch = forms.IntegerField(min_value=5, max_value=20, initial=10)


class PreprocessingForm(forms.Form):
    """预处理参数表单"""
    normMethod = forms.ChoiceField(choices=[
        ('quantile', 'Quantile Normalization'),
        ('swan', 'SWAN'),
        ('dasen', 'Dasen Normalization'),
        ('noob', 'NOOB'),
        ('none', 'None (raw beta values)')
    ], initial='quantile')
    batchMethod = forms.ChoiceField(choices=[
        ('combat', 'ComBat'),
        ('sva', 'SVA'),
        ('none', 'None')
    ], initial='combat')
    numSurrogateVars = forms.IntegerField(min_value=2, max_value=20, initial=5)
    valueTransformation = forms.ChoiceField(choices=[
        ('beta', 'Beta values (0-1)'),
        ('mvalue', 'M values'),
        ('both', 'Both beta and M values')
    ], initial='beta')
    probeMasking = forms.ChoiceField(choices=[
        ('cross-reactive', 'Exclude cross-reactive probes'),
        ('sex-chromosomes', 'Exclude sex chromosome probes'),
        ('both', 'Exclude both'),
        ('none', 'No probe masking')
    ], initial='both')
    smoothingMethod = forms.ChoiceField(choices=[
        ('none', 'None'),
        ('loess', 'LOESS smoothing'),
        ('gaussian', 'Gaussian kernel smoothing'),
        ('moving-average', 'Moving average')
    ], initial='none')
    smoothingWindow = forms.IntegerField(min_value=100, max_value=10000, initial=1000)


class DifferentialMethylationForm(forms.Form):
    """差异甲基化分析参数表单"""
    groupVariable = forms.CharField(initial='treatment')
    referenceGroup = forms.CharField(initial='control')
    testMethod = forms.ChoiceField(choices=[
        ('limma', 'Limma'),
        ('methylkit', 'MethylKit'),
        ('dss', 'DSS'),
        ('edgeR', 'edgeR')
    ], initial='limma')
    covariates = forms.CharField(required=False, initial='age,sex')
    pvalCutoff = forms.FloatField(min_value=0.001, max_value=0.1, step_size=0.01, initial=0.05)
    betaChangeCutoff = forms.FloatField(min_value=0.1, max_value=0.5, step_size=0.05, initial=0.2)
    mvalueChangeCutoff = forms.FloatField(min_value=0.2, max_value=2, step_size=0.1, initial=0.5)
    correctionMethod = forms.ChoiceField(choices=[
        ('fdr_bh', 'Benjamini-Hochberg (FDR)'),
        ('bonferroni', 'Bonferroni'),
        ('holm', 'Holm-Bonferroni'),
        ('qvalue', 'q-value')
    ], initial='fdr_bh')
    dmrMethod = forms.ChoiceField(choices=[
        ('bumphunter', 'Bumphunter'),
        ('dmrseq', 'DMRseq'),
        ('methylkit', 'MethylKit'),
        ('comb-p', 'Comb-p')
    ], initial='bumphunter')
    maxGap = forms.IntegerField(min_value=100, max_value=2000, initial=500)
    minCpGsDMR = forms.IntegerField(min_value=2, max_value=20, initial=3)
    minRegionLength = forms.IntegerField(min_value=10, max_value=1000, initial=50)


class RegionAnalysisForm(forms.Form):
    """区域分析参数表单"""
    predefinedRegions = forms.MultipleChoiceField(choices=[
        ('promoter', 'Promoters'),
        ('cgi', 'CpG Islands'),
        ('cgi-shore', 'CpG Island Shores'),
        ('cgi-shelf', 'CpG Island Shelves'),
        ('exon', 'Exons'),
        ('intron', 'Introns'),
        ('gene-body', 'Gene Bodies'),
        ('enhancer', 'Enhancers')
    ], initial=['promoter', 'cgi', 'cgi-shore'])
    customRegions = forms.ChoiceField(choices=[
        ('none', 'None'),
        ('uploaded', 'Use uploaded regions file')
    ], initial='none')
    promoterDefinition = forms.ChoiceField(choices=[
        ('tss-2kb-1kb', 'TSS ± 2kb'),
        ('tss-1kb-0kb', 'TSS -1kb to +0kb'),
        ('tss-0kb-1kb', 'TSS +0kb to +1kb'),
        ('tss-5kb-1kb', 'TSS -5kb to +1kb')
    ], initial='tss-2kb-1kb')
    geneAnnotation = forms.ChoiceField(choices=[
        ('ensembl', 'Ensembl'),
        ('refseq', 'RefSeq'),
        ('gencode', 'GENCODE')
    ], initial='ensembl')
    chromBinSize = forms.IntegerField(min_value=10, max_value=1000, initial=100)
    includeSexChromosomes = forms.ChoiceField(choices=[
        ('yes', 'Yes'),
        ('no', 'No')
    ], initial='yes')


class CorrelationForm(forms.Form):
    """相关性分析参数表单"""
    correlationMethod = forms.ChoiceField(choices=[
        ('pearson', 'Pearson correlation'),
        ('spearman', 'Spearman rank correlation'),
        ('both', 'Both')
    ], initial='pearson')
    correlationWindow = forms.IntegerField(min_value=1000, max_value=100000, initial=10000)
    minCorrelation = forms.FloatField(min_value=0.1, max_value=1.0, step_size=0.1, initial=0.3)
    correlationPval = forms.FloatField(min_value=0.001, max_value=0.1, step_size=0.01, initial=0.05)
    expressionType = forms.ChoiceField(choices=[
        ('mrna', 'mRNA expression'),
        ('transcript', 'Transcript expression'),
        ('protein', 'Protein levels')
    ], initial='mrna')
    expressionNorm = forms.ChoiceField(choices=[
        ('log2-tpm', 'Log2(TPM + 1)'),
        ('log2-fpkm', 'Log2(FPKM + 1)'),
        ('vst', 'Variance-stabilized'),
        ('none', 'None')
    ], initial='log2-tpm')
    simpleRegression = forms.BooleanField(initial=True, required=False)
    multipleRegression = forms.BooleanField(initial=False, required=False)
    segmentedRegression = forms.BooleanField(initial=False, required=False)


class EnrichmentForm(forms.Form):
    """富集分析参数表单"""
    geneSetSource = forms.ChoiceField(choices=[
        ('dmp-genes', 'Genes near DMPs'),
        ('dmr-genes', 'Genes overlapping DMRs'),
        ('correlated-genes', 'Methylation-expression correlated genes'),
        ('user-genes', 'User provided genes')
    ], initial='dmp-genes')
    geneDistance = forms.IntegerField(min_value=1, max_value=100, initial=10)
    enrichmentDb = forms.MultipleChoiceField(choices=[
        ('go-bp', 'Gene Ontology (Biological Process)'),
        ('go-cc', 'Gene Ontology (Cellular Component)'),
        ('go-mf', 'Gene Ontology (Molecular Function)'),
        ('kegg', 'KEGG Pathways'),
        ('reactome', 'Reactome Pathways'),
        ('msigdb', 'MSigDB Collections'),
        ('tfbs', 'TF Binding Sites')
    ], initial=['go-bp', 'go-cc', 'kegg'])
    enrichmentSpecies = forms.ChoiceField(choices=[
        ('human', 'Human'),
        ('mouse', 'Mouse'),
        ('rat', 'Rat'),
        ('zebrafish', 'Zebrafish')
    ], initial='human')
    enrichmentMethod = forms.ChoiceField(choices=[
        ('hypergeometric', 'Hypergeometric test'),
        ('fisher', 'Fisher\'s exact test'),
        ('gsea', 'GSEA')
    ], initial='hypergeometric')
    enrichPvalCutoff = forms.FloatField(min_value=0.001, max_value=0.1, step_size=0.01, initial=0.05)
    minGenesPerTerm = forms.IntegerField(min_value=2, max_value=50, initial=5)
    maxGenesPerTerm = forms.IntegerField(min_value=100, max_value=2000, initial=500)
