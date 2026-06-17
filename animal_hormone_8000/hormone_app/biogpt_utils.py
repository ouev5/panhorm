from transformers import BioGptTokenizer, BioGptForCausalLM
import torch

# 全局变量存储模型和分词器（避免重复加载）
MODEL = None
TOKENIZER = None

# 国内镜像地址（阿里云镜像，可根据需要替换为其他镜像）
HUGGINGFACE_MIRROR = "https://mirrors.aliyun.com/huggingface/models/"

def load_biogpt_model():
    """加载BioGPT模型和分词器（仅首次调用时加载）"""
    global MODEL, TOKENIZER
    if MODEL is None or TOKENIZER is None:
        # 加载预训练模型和分词器，指定国内镜像
        TOKENIZER = BioGptTokenizer.from_pretrained(
            "microsoft/biogpt",
            mirror=HUGGINGFACE_MIRROR  # 使用国内镜像
        )
        # 修复类名拼写错误（首字母大写）
        MODEL = BioGptForCausalLM.from_pretrained(
            "microsoft/biogpt",
            mirror=HUGGINGFACE_MIRROR  # 使用国内镜像
        )
        
        # 如果有GPU，将模型移至GPU
        if torch.cuda.is_available():
            MODEL = MODEL.to("cuda")
            MODEL.eval()  # 推理模式，关闭dropout等训练相关层
    return MODEL, TOKENIZER

def generate_biomed_response(prompt, max_length=200):
    """
    调用BioGPT生成生物医学相关回答
    prompt: 输入的问题或提示文本（如"Describe the function of BRCA1 gene"）
    max_length: 生成文本的最大长度
    """
    model, tokenizer = load_biogpt_model()
    
    # 对输入文本进行编码
    inputs = tokenizer(prompt, return_tensors="pt", truncation=True, padding=True)
    
    # 如果有GPU，将输入移至GPU
    if torch.cuda.is_available():
        inputs = {k: v.to("cuda") for k, v in inputs.items()}
    
    # 生成回答（设置no_repeat_ngram_size避免重复）
    with torch.no_grad():  # 关闭梯度计算，节省内存
        outputs = model.generate(
            **inputs,
            max_length=max_length,
            no_repeat_ngram_size=2,
            do_sample=True,
            temperature=0.7,  # 控制生成多样性（值越低越确定）
            pad_token_id=tokenizer.eos_token_id  # 修复可能的padding警告
        )
    
    # 解码生成的文本
    response = tokenizer.decode(outputs[0], skip_special_tokens=True)
    return response
