from django.core.management.base import BaseCommand
# 正确导入路径：hormone_app.util.rag_retriever（文件在util目录下）
from hormone_app.utils.rag_retriever import vector_store

class Command(BaseCommand):
    help = '手动构建激素数据语义检索向量索引'

    def handle(self, *args, **options):
        self.stdout.write(self.style.SUCCESS('开始手动构建语义检索向量索引...'))
        # 临时关闭强制禁止开关
        vector_store.FORCE_DISABLE_AUTO_BUILD = False
        # 执行构建
        vector_store.load_data_to_vector()
        # 重新打开强制开关
        vector_store.FORCE_DISABLE_AUTO_BUILD = True
        # 构建结果提示
        if vector_store.is_loaded:
            self.stdout.write(self.style.SUCCESS(f'✅ 向量索引构建完成！共包含 {len(vector_store.texts)} 条记录'))
        else:
            self.stdout.write(self.style.ERROR('❌ 向量索引构建失败'))
