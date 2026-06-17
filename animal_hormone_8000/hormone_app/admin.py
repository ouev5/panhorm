from django.contrib import admin
from .models import (
    Hormone, 
    Literature, 
    HormoneRelatedGene, 
    HormoneReceptorInfo, 
    HormoneReceptorFull
)


@admin.register(Hormone)
class HormoneAdmin(admin.ModelAdmin):
    # 结合数据导入和模型，假设实际有这些字段（需和Hormone模型实际字段一致，这里按导入逻辑推测常用展示字段）
    list_display = (
        'id',
        'receptor_name',
        'gene_family',
        'gene_name_symbol',
        'ligand',
        'ligand_type_comments',
    )
    search_fields = ('receptor_name', 'gene_name_symbol', 'ligand')
    list_filter = ('gene_family',)  # 若有receptor_type实际字段再加上，这里先按合理推测
    ordering = ('receptor_name',)


@admin.register(Literature)
class LiteratureAdmin(admin.ModelAdmin):
    list_display = ('title', 'authors', 'year', 'lit_type')
    list_filter = ('lit_type', 'year')
    search_fields = ('title', 'authors')  # 若模型无abstract字段，去掉该搜索项
    ordering = ('-year', 'title')  
    #date_hierarchy = 'year'  


@admin.register(HormoneRelatedGene)
class HormoneRelatedGeneAdmin(admin.ModelAdmin):
    # 对应HormoneRelatedGene模型实际字段，从导入逻辑看有这些核心字段
    list_display = ('id', 'hormone_name', 'related_genes', 'pmid', 'gene_sequence')  
    search_fields = ('hormone_name', 'related_genes', 'gene_sequence')  
    list_filter = ('hormone_name',)  # 可按实际需求加更多筛选，比如pmid范围等
    # 按名称排序，方便查找
    ordering = ('hormone_name',)  


@admin.register(HormoneReceptorInfo)
class HormoneReceptorInfoAdmin(admin.ModelAdmin):
    # 匹配HormoneReceptorInfo模型，用导入涉及的实际字段
    list_display = ('id', 'hormone_name', 'receptor_uniprot_id', 'receptor_name', 'receptor_species_name')  
    search_fields = ('hormone_name', 'receptor_uniprot_id', 'receptor_name')  
    list_filter = ('receptor_species_name',)  


@admin.register(HormoneReceptorFull)
class HormoneReceptorFullAdmin(admin.ModelAdmin):
    # 对应HormoneReceptorFull模型导入的字段，展示关键信息
    list_display = ('id', 'hormone_name', 'hormone_species_name', 'receptor_name', 'receptor_species_name')  
    search_fields = ('hormone_name', 'receptor_name', 'hormone_coding_genes_sequence')  
    list_filter = ('hormone_species_name', 'receptor_species_name')  