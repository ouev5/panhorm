import json
import hashlib
import time
from typing import List, Dict, Any, Optional
from models import UserRequest, PipelineTemplate, PipelineStep
from rag_retriever import RAGRetriever
from llm_client import LLMPlannerClient
from registry import PipelineRegistry


class Planner:
    def __init__(self, config: dict):
        """初始化规划器
        
        Args:
            config: 配置字典
        """
        self.config = config
        self.rag_retriever = RAGRetriever(config.get('rag_persist_dir', './rag_index'))
        self.llm_client = LLMPlannerClient(config.get('llm', {}))
        self.registry = PipelineRegistry()
        self.cache = {}  # 计划缓存
        self.cache_ttl = config.get('cache_ttl', 3600)  # 缓存过期时间（秒）
        self.max_retries = config.get('max_retries', 3)  # 最大重试次数

    def _get_request_hash(self, user_request: UserRequest) -> str:
        """计算用户请求的哈希值，用于缓存键
        
        Args:
            user_request: 用户请求
            
        Returns:
            str: 哈希值
        """
        request_data = {
            'data_paths': user_request.data_paths,
            'data_description': user_request.data_description,
            'goal': user_request.goal,
            'task_type': user_request.task_type.value
        }
        return hashlib.md5(json.dumps(request_data, sort_keys=True).encode()).hexdigest()

    def _get_from_cache(self, user_request: UserRequest) -> Optional[Dict[str, Any]]:
        """从缓存中获取计划
        
        Args:
            user_request: 用户请求
            
        Returns:
            Optional[Dict[str, Any]]: 缓存的计划数据，如果不存在则返回None
        """
        cache_key = self._get_request_hash(user_request)
        if cache_key in self.cache:
            cached_data = self.cache[cache_key]
            if time.time() - cached_data['timestamp'] < self.cache_ttl:
                return cached_data
        return None

    def _save_to_cache(self, user_request: UserRequest, plan_data: Dict[str, Any]):
        """将计划保存到缓存
        
        Args:
            user_request: 用户请求
            plan_data: 计划数据
        """
        cache_key = self._get_request_hash(user_request)
        self.cache[cache_key] = {
            'timestamp': time.time(),
            'plan': plan_data
        }

    def _retrieve_context(self, user_request: UserRequest, top_k: int = 5) -> List[Dict[str, Any]]:
        """检索相关文档
        
        Args:
            user_request: 用户请求
            top_k: 返回的文档数量
            
        Returns:
            List[Dict[str, Any]]: 检索到的文档
        """
        # 构建查询语句
        query = f"{user_request.goal} {user_request.data_description}"
        
        # 检索相关文档
        try:
            documents = self.rag_retriever.retrieve(query, top_k=top_k)
            return [{
                'content': doc.page_content,
                'source': doc.metadata.get('source', '')
            } for doc in documents]
        except Exception as e:
            print(f"RAG retrieval error: {e}")
            return []

    def _build_prompt(self, user_request: UserRequest, context: List[Dict[str, Any]]) -> List[Dict[str, str]]:
        """构建Prompt
        
        Args:
            user_request: 用户请求
            context: 检索到的上下文
            
        Returns:
            List[Dict[str, str]]: Prompt消息列表
        """
        # 获取可用工具列表
        available_tools = [
            "Trimmomatic", "HISAT2", "samtools", "featureCounts", "HTSeq", "DESeq2",
            "Scanpy", "MACS2", "GATK", "MaxQuant", "OpenMS"
        ]
        
        # 构建上下文字符串
        context_str = "\n".join([f"[Source: {doc['source']}]\n{doc['content'][:500]}..." for doc in context])
        
        # 构建系统消息
        system_message = {
            "role": "system",
            "content": "你是一个生物信息学分析专家，负责根据用户请求生成详细的分析计划。请根据用户提供的信息、检索到的上下文和可用工具，生成一个完整的分析流程计划。"
        }
        
        # 构建用户消息
        user_message = {
            "role": "user",
            "content": f"用户请求：\n"\
                      f"数据路径：{user_request.data_paths}\n"\
                      f"数据描述：{user_request.data_description}\n"\
                      f"分析目标：{user_request.goal}\n"\
                      f"任务类型：{user_request.task_type.value}\n"\
                      f"\n可用工具：{', '.join(available_tools)}\n"\
                      f"\n检索到的相关文档：\n{context_str}\n"\
                      f"\n请生成一个JSON格式的分析计划，包含以下字段：\n"\
                      f"- name: 计划名称\n"\
                      f"- description: 计划描述\n"\
                      f"- task_type: 任务类型（deterministic或ambiguous）\n"\
                      f"- required_inputs: 必需的输入文件列表\n"\
                      f"- output_patterns: 输出文件模式列表\n"\
                      f"- steps: 步骤列表，每个步骤包含：\n"\
                      f"  - step_id: 步骤ID（从1开始）\n"\
                      f"  - tool_name: 工具名称\n"\
                      f"  - tool_version: 工具版本\n"\
                      f"  - command_template: 命令模板，使用{变量}占位符\n"\
                      f"  - input_files: 输入文件列表\n"\
                      f"  - output_files: 输出文件列表\n"\
                      f"  - parameters: 参数字典\n"
        }
        
        return [system_message, user_message]

    def _parse_llm_response(self, response: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """解析LLM响应
        
        Args:
            response: LLM响应
            
        Returns:
            Optional[Dict[str, Any]]: 解析后的计划数据
        """
        try:
            # 提取内容
            content = response['choices'][0]['message']['content']
            
            # 提取JSON部分
            import re
            json_match = re.search(r'\{[\s\S]*\}', content)
            if not json_match:
                return None
            
            json_str = json_match.group(0)
            plan_data = json.loads(json_str)
            
            # 验证必要字段
            required_fields = ['name', 'description', 'task_type', 'required_inputs', 'output_patterns', 'steps']
            for field in required_fields:
                if field not in plan_data:
                    return None
            
            return plan_data
        except Exception as e:
            print(f"Failed to parse LLM response: {e}")
            return None

    def _validate_plan(self, plan_data: Dict[str, Any]) -> bool:
        """验证计划
        
        Args:
            plan_data: 计划数据
            
        Returns:
            bool: 验证是否通过
        """
        try:
            # 验证步骤ID是否连续
            step_ids = [step['step_id'] for step in plan_data['steps']]
            if sorted(step_ids) != list(range(1, len(step_ids) + 1)):
                return False
            
            # 验证每个步骤的必需字段
            for step in plan_data['steps']:
                required_step_fields = ['step_id', 'tool_name', 'tool_version', 'command_template', 'input_files', 'output_files', 'parameters']
                for field in required_step_fields:
                    if field not in step:
                        return False
            
            return True
        except Exception:
            return False

    def _convert_to_pipeline_template(self, plan_data: Dict[str, Any]) -> PipelineTemplate:
        """将计划数据转换为PipelineTemplate
        
        Args:
            plan_data: 计划数据
            
        Returns:
            PipelineTemplate: 管道模板
        """
        # 转换步骤
        steps = []
        for step_data in plan_data['steps']:
            step = PipelineStep(
                step_id=step_data['step_id'],
                tool_name=step_data['tool_name'],
                tool_version=step_data['tool_version'],
                command_template=step_data['command_template'],
                input_files=step_data['input_files'],
                output_files=step_data['output_files'],
                parameters=step_data['parameters']
            )
            steps.append(step)
        
        # 创建PipelineTemplate
        template = PipelineTemplate(
            name=plan_data['name'],
            description=plan_data['description'],
            task_type=plan_data['task_type'],
            steps=steps,
            required_inputs=plan_data['required_inputs'],
            output_patterns=plan_data['output_patterns']
        )
        
        return template

    async def plan(self, user_request: UserRequest) -> PipelineTemplate:
        """生成分析计划
        
        Args:
            user_request: 用户请求
            
        Returns:
            PipelineTemplate: 分析计划
        """
        # 检查缓存
        cached_data = self._get_from_cache(user_request)
        if cached_data:
            return self._convert_to_pipeline_template(cached_data['plan'])
        
        # 多轮尝试
        for retry in range(self.max_retries):
            # 检索上下文
            context = self._retrieve_context(user_request)
            
            # 构建Prompt
            messages = self._build_prompt(user_request, context)
            
            # 调用LLM
            response = await self.llm_client.generate(messages)
            
            # 解析响应
            plan_data = self._parse_llm_response(response)
            
            # 验证计划
            if plan_data and self._validate_plan(plan_data):
                # 保存到缓存
                self._save_to_cache(user_request, plan_data)
                
                # 转换为PipelineTemplate
                return self._convert_to_pipeline_template(plan_data)
            
            print(f"Plan validation failed, retrying ({retry + 1}/{self.max_retries})...")
        
        # 回退到规则
        return self._fallback_plan(user_request)

    def _fallback_plan(self, user_request: UserRequest) -> PipelineTemplate:
        """回退计划（当LLM生成失败时使用）
        
        Args:
            user_request: 用户请求
            
        Returns:
            PipelineTemplate: 回退计划
        """
        # 根据任务类型生成回退计划
        if 'RNA-seq' in user_request.data_description or '转录组' in user_request.data_description:
            # 回退到RNA-seq差异表达流程
            return self.registry.get('rnaseq_diffexpr')
        elif '单细胞' in user_request.data_description:
            # 回退到单细胞聚类流程
            return self.registry.get('scrna_cluster')
        elif 'ChIP-seq' in user_request.data_description:
            # 回退到ChIP-seq峰值检测流程
            return self.registry.get('chipseq_peaks')
        else:
            # 默认回退计划
            steps = [
                PipelineStep(
                    step_id=1,
                    tool_name='echo',
                    tool_version='1.0',
                    command_template='echo "Default plan executed" > {output_dir}/result.txt',
                    input_files=[],
                    output_files=['{output_dir}/result.txt'],
                    parameters={}
                )
            ]
            
            return PipelineTemplate(
                name='default_plan',
                description='默认回退计划',
                task_type='deterministic',
                steps=steps,
                required_inputs=[],
                output_patterns=['{output_dir}/result.txt']
            )
