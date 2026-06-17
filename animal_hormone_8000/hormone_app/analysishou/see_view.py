# /www/wwwroot/default/animal_hormone/hormone_app/analysishou/see_view.py

import os
import uuid
import requests
import base64
from django.conf import settings
from django.http import JsonResponse
from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_POST
from datetime import datetime
from pdf2image import convert_from_path
from PIL import Image
import io
import tempfile

# 配置文件上传路径
PAPER_UPLOAD_DIR = os.path.join(settings.MEDIA_ROOT, 'papers')
os.makedirs(PAPER_UPLOAD_DIR, exist_ok=True)

def page11(request):
    """See页面 - 论文内容提取"""
    return render(request, 'page1.html')

@require_POST
@csrf_exempt
def paper_qa(request):
    """
    处理论文问答接口
    接收PDF/TXT/图片文件，转换为千问视觉模型支持的格式
    """
    print("="*50)
    print("收到 paper_qa 请求")
    print("="*50)
    
    try:
        # 获取上传的文件和问题
        uploaded_file = request.FILES.get('file')
        question = request.POST.get('question')
        model_name = request.POST.get('model', 'qwen-vl-max')
        
        if not uploaded_file:
            return JsonResponse({
                'status': 'error',
                'message': '没有收到文件'
            }, status=400)
            
        if not question:
            return JsonResponse({
                'status': 'error',
                'message': '没有收到问题'
            }, status=400)
        
        # 验证文件类型
        file_ext = os.path.splitext(uploaded_file.name)[1].lower()
        allowed_extensions = ['.pdf', '.txt', '.jpg', '.jpeg', '.png', '.bmp']
        if file_ext not in allowed_extensions:
            return JsonResponse({
                'status': 'error',
                'message': f'不支持的文件类型: {file_ext}，支持的类型: PDF, TXT, 图片'
            }, status=400)
        
        # 生成唯一文件名并保存
        file_id = str(uuid.uuid4())[:8]
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        safe_filename = f"{timestamp}_{file_id}_{uploaded_file.name}"
        file_path = os.path.join(PAPER_UPLOAD_DIR, safe_filename)
        
        # 保存文件
        print(f"保存文件到: {file_path}")
        with open(file_path, 'wb+') as destination:
            for chunk in uploaded_file.chunks():
                destination.write(chunk)
        
        # 根据文件类型处理
        if file_ext == '.pdf':
            # PDF文件：转换为图片再调用API
            answer = process_pdf_and_call_api(file_path, question, model_name)
        elif file_ext in ['.txt']:
            # TXT文件：直接作为文本处理
            answer = process_txt_and_call_api(file_path, question, model_name)
        else:
            # 图片文件：直接调用API
            answer = call_qwen_vl_with_image(file_path, question, model_name)
        
        return JsonResponse({
            'status': 'success',
            'answer': answer,
            'file_name': uploaded_file.name
        })
        
    except Exception as e:
        print(f"处理异常: {str(e)}")
        import traceback
        traceback.print_exc()
        return JsonResponse({
            'status': 'error',
            'message': f'处理失败: {str(e)}'
        }, status=500)

def pdf_to_images(pdf_path, max_pages=5):
    """
    将PDF转换为图片，限制页数避免token超限
    """
    try:
        print(f"转换PDF为图片: {pdf_path}")
        images = convert_from_path(pdf_path, dpi=150, fmt='jpeg')
        
        # 限制页数，避免API调用超时
        if len(images) > max_pages:
            print(f"PDF页数过多 ({len(images)}页)，只处理前{max_pages}页")
            images = images[:max_pages]
        
        # 将图片转换为base64
        image_data_list = []
        for i, image in enumerate(images):
            # 压缩图片大小
            img_buffer = io.BytesIO()
            image.save(img_buffer, format='JPEG', quality=85, optimize=True)
            img_base64 = base64.b64encode(img_buffer.getvalue()).decode('utf-8')
            image_data_list.append({
                'index': i,
                'data': img_base64,
                'mime': 'image/jpeg'
            })
            print(f"第{i+1}页转换完成，大小: {len(img_base64)/1024:.1f}KB")
        
        return image_data_list
    except Exception as e:
        print(f"PDF转换失败: {str(e)}")
        raise Exception(f"PDF转换失败: {str(e)}")

def process_pdf_and_call_api(pdf_path, question, model_name):
    """
    处理PDF文件：转换为图片后调用API
    """
    api_key = getattr(settings, 'DASHSCOPE_API_KEY', '')
    
    if not api_key or api_key == 'YOUR_DASHSCOPE_API_KEY':
        return f"""【模拟回答】\n\n您上传的PDF: {os.path.basename(pdf_path)}\n问题: {question}\n\n已转换为图片格式，但由于未配置API密钥，这是模拟回答。"""
    
    try:
        # 将PDF转换为图片
        image_list = pdf_to_images(pdf_path, max_pages=3)  # 限制3页，避免token过多
        
        # 构建多模态消息
        import json
        
        headers = {
            'Authorization': f'Bearer {api_key}',
            'Content-Type': 'application/json'
        }
        
        # 构建content数组，包含所有图片和问题
        content = []
        
        # 添加图片（最多3张）
        for img in image_list:
            content.append({
                'image': f'data:{img["mime"]};base64,{img["data"]}'
            })
        
        # 添加文本问题
        content.append({
            'text': f'这是PDF文档的第1-{len(image_list)}页内容。请基于这些页面回答以下问题：{question}\n\n请详细分析，如果信息不足，请说明。'
        })
        
        messages = [
            {
                'role': 'user',
                'content': content
            }
        ]
        
        payload = {
            'model': 'qwen-vl-max',  # 或 'qwen-vl-plus'
            'input': {
                'messages': messages
            },
            'parameters': {
                'temperature': 0.7,
                'max_tokens': 2000
            }
        }
        
        print("调用DashScope API处理PDF转换后的图片...")
        
        response = requests.post(
            'https://dashscope.aliyuncs.com/api/v1/services/aigc/multimodal-generation/generation',
            headers=headers,
            json=payload,
            timeout=180  # 增加超时时间
        )
        
        print(f"API响应状态: {response.status_code}")
        
        if response.status_code == 200:
            result = response.json()
            if 'output' in result and 'choices' in result['output']:
                return result['output']['choices'][0]['message']['content']
            elif 'output' in result and 'text' in result['output']:
                return result['output']['text']
            else:
                return f"API返回格式异常: {json.dumps(result, ensure_ascii=False)[:500]}"
        else:
            error_text = response.text
            print(f"API错误: {error_text}")
            return f"API调用失败 (HTTP {response.status_code}): {error_text[:500]}"
            
    except Exception as e:
        print(f"PDF处理异常: {str(e)}")
        import traceback
        traceback.print_exc()
        return f"PDF处理失败: {str(e)}\n\n已保存PDF文件，但无法调用API。您可以稍后重试。"

def process_txt_and_call_api(txt_path, question, model_name):
    """
    处理TXT文件：直接作为文本调用API
    """
    api_key = getattr(settings, 'DASHSCOPE_API_KEY', '')
    
    try:
        # 读取TXT文件内容
        with open(txt_path, 'r', encoding='utf-8', errors='ignore') as f:
            text_content = f.read()[:8000]  # 限制长度
        
        headers = {
            'Authorization': f'Bearer {api_key}',
            'Content-Type': 'application/json'
        }
        
        messages = [
            {
                'role': 'user',
                'content': f'以下是文档内容：\n\n{text_content}\n\n请基于以上内容回答问题：{question}\n\n请给出详细、准确的回答。'
            }
        ]
        
        payload = {
            'model': 'qwen-max',  # 对纯文本可以使用qwen-max
            'input': {
                'messages': messages
            },
            'parameters': {
                'temperature': 0.7,
                'max_tokens': 2000
            }
        }
        
        response = requests.post(
            'https://dashscope.aliyuncs.com/api/v1/services/aigc/text-generation/generation',
            headers=headers,
            json=payload,
            timeout=120
        )
        
        if response.status_code == 200:
            result = response.json()
            if 'output' in result and 'choices' in result['output']:
                return result['output']['choices'][0]['message']['content']
            elif 'output' in result and 'text' in result['output']:
                return result['output']['text']
        
        return f"TXT处理失败: {response.text}"
        
    except Exception as e:
        return f"TXT处理异常: {str(e)}"

def call_qwen_vl_with_image(image_path, question, model_name):
    """
    直接处理图片文件
    """
    api_key = getattr(settings, 'DASHSCOPE_API_KEY', '')
    
    try:
        with open(image_path, 'rb') as f:
            img_base64 = base64.b64encode(f.read()).decode('utf-8')
        
        headers = {
            'Authorization': f'Bearer {api_key}',
            'Content-Type': 'application/json'
        }
        
        messages = [
            {
                'role': 'user',
                'content': [
                    {
                        'image': f'data:image/jpeg;base64,{img_base64}'
                    },
                    {
                        'text': f'请分析这张图片并回答：{question}'
                    }
                ]
            }
        ]
        
        payload = {
            'model': model_name,
            'input': {
                'messages': messages
            },
            'parameters': {
                'temperature': 0.7,
                'max_tokens': 2000
            }
        }
        
        response = requests.post(
            'https://dashscope.aliyuncs.com/api/v1/services/aigc/multimodal-generation/generation',
            headers=headers,
            json=payload,
            timeout=120
        )
        
        if response.status_code == 200:
            result = response.json()
            if 'output' in result and 'choices' in result['output']:
                return result['output']['choices'][0]['message']['content']
        
        return f"图片处理失败: {response.text}"
        
    except Exception as e:
        return f"图片处理异常: {str(e)}"