import asyncio
import time
import json
from typing import List, Dict, Any, Optional, AsyncGenerator


class LLMPlannerClient:
    def __init__(self, config: dict):
        """初始化LLM规划客户端
        
        Args:
            config: 配置字典
        """
        self.config = config
        self.mode = config.get('mode', 'mock')  # online, local, mock
        self.max_retries = config.get('max_retries', 3)
        self.base_delay = config.get('base_delay', 1)  # 基础重试延迟（秒）
        
        # Token计数和成本跟踪
        self.total_tokens = 0
        self.total_cost = 0.0
        
        # 模型配置
        if self.mode == 'online':
            self.online_config = config.get('online', {})
            self.model = self.online_config.get('model', 'gpt-4-turbo-preview')
            self.temperature = self.online_config.get('temperature', 0.2)
        elif self.mode == 'local':
            self.local_config = config.get('local', {})
            self.api_url = self.local_config.get('api_url', 'http://localhost:8000')
            self.model_name = self.local_config.get('model_name', 'deepseek-coder-6.7b')
        
        # 成本配置（美元/1000 tokens）
        self.cost_config = {
            'gpt-4': {'input': 0.03, 'output': 0.06},
            'gpt-4-turbo-preview': {'input': 0.01, 'output': 0.03},
            'gpt-3.5-turbo': {'input': 0.0015, 'output': 0.002},
        }

    async def generate(self, messages: List[Dict[str, str]], **kwargs) -> Dict[str, Any]:
        """生成响应
        
        Args:
            messages: 消息列表
            **kwargs: 额外参数
            
        Returns:
            Dict[str, Any]: 生成结果
        """
        retry_count = 0
        while retry_count < self.max_retries:
            try:
                if self.mode == 'online':
                    return await self._generate_online(messages, **kwargs)
                elif self.mode == 'local':
                    return await self._generate_local(messages, **kwargs)
                else:  # mock
                    return await self._generate_mock(messages, **kwargs)
            except Exception as e:
                retry_count += 1
                if retry_count >= self.max_retries:
                    raise
                
                # 指数退避
                delay = self.base_delay * (2 ** (retry_count - 1))
                print(f"Error: {e}. Retrying in {delay}s...")
                await asyncio.sleep(delay)

    async def _generate_online(self, messages: List[Dict[str, str]], **kwargs) -> Dict[str, Any]:
        """使用OpenAI API生成响应
        
        Args:
            messages: 消息列表
            **kwargs: 额外参数
            
        Returns:
            Dict[str, Any]: 生成结果
        """
        import openai
        
        # 获取API密钥
        api_key = self.online_config.get('api_key') or kwargs.get('api_key')
        if not api_key:
            raise ValueError("OpenAI API key is required")
        
        openai.api_key = api_key
        
        # 调用API
        response = await openai.ChatCompletion.acreate(
            model=self.model,
            messages=messages,
            temperature=self.temperature,
            **kwargs
        )
        
        # 跟踪Token和成本
        usage = response.get('usage', {})
        input_tokens = usage.get('prompt_tokens', 0)
        output_tokens = usage.get('completion_tokens', 0)
        self._track_usage(input_tokens, output_tokens)
        
        return response

    async def _generate_local(self, messages: List[Dict[str, str]], **kwargs) -> Dict[str, Any]:
        """使用本地LLM服务生成响应
        
        Args:
            messages: 消息列表
            **kwargs: 额外参数
            
        Returns:
            Dict[str, Any]: 生成结果
        """
        import httpx
        
        # 构建请求
        payload = {
            'model': self.model_name,
            'messages': messages,
            'temperature': kwargs.get('temperature', self.temperature),
            **kwargs
        }
        
        # 发送请求
        async with httpx.AsyncClient() as client:
            response = await client.post(
                f"{self.api_url}/v1/chat/completions",
                json=payload
            )
            response.raise_for_status()
            result = response.json()
        
        # 跟踪Token（如果返回了使用情况）
        usage = result.get('usage', {})
        input_tokens = usage.get('prompt_tokens', 0)
        output_tokens = usage.get('completion_tokens', 0)
        self._track_usage(input_tokens, output_tokens)
        
        return result

    async def _generate_mock(self, messages: List[Dict[str, str]], **kwargs) -> Dict[str, Any]:
        """生成模拟响应
        
        Args:
            messages: 消息列表
            **kwargs: 额外参数
            
        Returns:
            Dict[str, Any]: 模拟结果
        """
        # 模拟响应
        mock_response = {
            'id': 'mock-response-123',
            'object': 'chat.completion',
            'created': int(time.time()),
            'model': 'mock-model',
            'choices': [{
                'index': 0,
                'message': {
                    'role': 'assistant',
                    'content': '这是一个模拟响应，用于测试。'
                },
                'finish_reason': 'stop'
            }],
            'usage': {
                'prompt_tokens': 10,
                'completion_tokens': 20,
                'total_tokens': 30
            }
        }
        
        # 跟踪Token
        self._track_usage(10, 20)
        
        # 模拟延迟
        await asyncio.sleep(0.5)
        
        return mock_response

    async def stream(self, messages: List[Dict[str, str]], **kwargs) -> AsyncGenerator[Dict[str, Any], None]:
        """流式生成响应
        
        Args:
            messages: 消息列表
            **kwargs: 额外参数
            
        Yields:
            Dict[str, Any]: 流式响应数据
        """
        if self.mode == 'online':
            async for chunk in self._stream_online(messages, **kwargs):
                yield chunk
        elif self.mode == 'local':
            async for chunk in self._stream_local(messages, **kwargs):
                yield chunk
        else:  # mock
            async for chunk in self._stream_mock(messages, **kwargs):
                yield chunk

    async def _stream_online(self, messages: List[Dict[str, str]], **kwargs) -> AsyncGenerator[Dict[str, Any], None]:
        """使用OpenAI API流式生成响应
        
        Args:
            messages: 消息列表
            **kwargs: 额外参数
            
        Yields:
            Dict[str, Any]: 流式响应数据
        """
        import openai
        
        # 获取API密钥
        api_key = self.online_config.get('api_key') or kwargs.get('api_key')
        if not api_key:
            raise ValueError("OpenAI API key is required")
        
        openai.api_key = api_key
        
        # 调用API
        async for chunk in await openai.ChatCompletion.acreate(
            model=self.model,
            messages=messages,
            temperature=self.temperature,
            stream=True,
            **kwargs
        ):
            yield chunk

    async def _stream_local(self, messages: List[Dict[str, str]], **kwargs) -> AsyncGenerator[Dict[str, Any], None]:
        """使用本地LLM服务流式生成响应
        
        Args:
            messages: 消息列表
            **kwargs: 额外参数
            
        Yields:
            Dict[str, Any]: 流式响应数据
        """
        import httpx
        
        # 构建请求
        payload = {
            'model': self.model_name,
            'messages': messages,
            'temperature': kwargs.get('temperature', self.temperature),
            'stream': True,
            **kwargs
        }
        
        # 发送请求
        async with httpx.AsyncClient(timeout=None) as client:
            async with client.stream(
                'POST',
                f"{self.api_url}/v1/chat/completions",
                json=payload
            ) as response:
                response.raise_for_status()
                async for chunk in response.aiter_text():
                    if chunk.strip():
                        # 处理SSE格式
                        for line in chunk.split('\n'):
                            line = line.strip()
                            if line.startswith('data: '):
                                data = line[6:]
                                if data == '[DONE]':
                                    yield {'done': True}
                                else:
                                    try:
                                        yield json.loads(data)
                                    except json.JSONDecodeError:
                                        pass

    async def _stream_mock(self, messages: List[Dict[str, str]], **kwargs) -> AsyncGenerator[Dict[str, Any], None]:
        """生成模拟流式响应
        
        Args:
            messages: 消息列表
            **kwargs: 额外参数
            
        Yields:
            Dict[str, Any]: 模拟流式响应数据
        """
        # 模拟流式响应
        mock_chunks = [
            {
                'id': 'mock-stream-123',
                'object': 'chat.completion.chunk',
                'created': int(time.time()),
                'model': 'mock-model',
                'choices': [{
                    'index': 0,
                    'delta': {'role': 'assistant'},
                    'finish_reason': None
                }]
            },
            {
                'id': 'mock-stream-123',
                'object': 'chat.completion.chunk',
                'created': int(time.time()),
                'model': 'mock-model',
                'choices': [{
                    'index': 0,
                    'delta': {'content': '这是一个模拟的流式响应，'},
                    'finish_reason': None
                }]
            },
            {
                'id': 'mock-stream-123',
                'object': 'chat.completion.chunk',
                'created': int(time.time()),
                'model': 'mock-model',
                'choices': [{
                    'index': 0,
                    'delta': {'content': '用于测试流式生成功能。'},
                    'finish_reason': None
                }]
            },
            {
                'id': 'mock-stream-123',
                'object': 'chat.completion.chunk',
                'created': int(time.time()),
                'model': 'mock-model',
                'choices': [{
                    'index': 0,
                    'delta': {},
                    'finish_reason': 'stop'
                }]
            }
        ]
        
        # 模拟流式输出
        for chunk in mock_chunks:
            await asyncio.sleep(0.2)
            yield chunk

    def _track_usage(self, input_tokens: int, output_tokens: int):
        """跟踪Token使用和成本
        
        Args:
            input_tokens: 输入Token数
            output_tokens: 输出Token数
        """
        self.total_tokens += input_tokens + output_tokens
        
        # 计算成本
        model_cost = self.cost_config.get(self.model, {'input': 0, 'output': 0})
        input_cost = (input_tokens / 1000) * model_cost['input']
        output_cost = (output_tokens / 1000) * model_cost['output']
        self.total_cost += input_cost + output_cost

    def get_usage_stats(self) -> Dict[str, Any]:
        """获取使用统计信息
        
        Returns:
            Dict[str, Any]: 使用统计信息
        """
        return {
            'total_tokens': self.total_tokens,
            'total_cost': round(self.total_cost, 4),
            'mode': self.mode
        }
