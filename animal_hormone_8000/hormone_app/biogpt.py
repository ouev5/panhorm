from django.apps import AppConfig
import torch
from transformers import AutoTokenizer, AutoModelForCausalLM

class HormoneAppConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'hormone_app'
    
    # 模型全局实例
    tokenizer = None
    model = None

    def ready(self):
        # 模型路径（替换为你的实际路径）
        model_path = "/root/BioGPT/bioGPT_model"
        
        # 仅在项目启动时加载一次
        try:
            self.tokenizer = AutoTokenizer.from_pretrained(model_path)
            self.model = AutoModelForCausalLM.from_pretrained(model_path)
            
            # 启用GPU加速（如果有）
            if torch.cuda.is_available():
                self.model = self.model.to("cuda")
        except Exception as e:
            print(f"模型加载失败: {str(e)}")