import json
import hashlib
import time
import asyncio
import os
import pandas as pd
import numpy as np
from typing import Optional, Dict, Any, List
from models import ValidationResult


class LLMValidator:
    def __init__(self, config: dict):
        """初始化LLM验证器
        
        Args:
            config: 配置字典，包含模型选择、API密钥等
        """
        self.config = config
        self.model_type = config.get('model_type', 'openai')  # 'openai' or 'local'
        self.openai_api_key = config.get('openai_api_key')
        self.local_model_url = config.get('local_model_url', 'http://localhost:8000')
        self.timeout = config.get('timeout', 30)
        self.cache = {}
        self.cache_ttl = config.get('cache_ttl', 3600)  # 缓存过期时间（秒）

    def _get_content_hash(self, content: str) -> str:
        """计算内容的哈希值，用于缓存键
        
        Args:
            content: 要哈希的内容
            
        Returns:
            str: 哈希值
        """
        return hashlib.md5(content.encode()).hexdigest()

    def _get_from_cache(self, content: str, task_name: str) -> Optional[ValidationResult]:
        """从缓存中获取验证结果
        
        Args:
            content: 要验证的内容
            task_name: 任务名称
            
        Returns:
            Optional[ValidationResult]: 缓存的验证结果，如果不存在则返回None
        """
        cache_key = f"{task_name}:{self._get_content_hash(content)}"
        if cache_key in self.cache:
            cached_data = self.cache[cache_key]
            if time.time() - cached_data['timestamp'] < self.cache_ttl:
                return ValidationResult(**cached_data['result'])
            else:
                # 缓存过期，删除
                del self.cache[cache_key]
        return None

    def _save_to_cache(self, content: str, task_name: str, result: ValidationResult):
        """将验证结果保存到缓存
        
        Args:
            content: 要验证的内容
            task_name: 任务名称
            result: 验证结果
        """
        cache_key = f"{task_name}:{self._get_content_hash(content)}"
        self.cache[cache_key] = {
            'timestamp': time.time(),
            'result': result.to_dict()
        }

    async def _call_openai(self, prompt: str) -> Optional[str]:
        """调用OpenAI API
        
        Args:
            prompt: 提示词
            
        Returns:
            Optional[str]: API响应
        """
        try:
            import openai
            openai.api_key = self.openai_api_key
            
            async with asyncio.timeout(self.timeout):
                response = await openai.ChatCompletion.acreate(
                    model="gpt-4",
                    messages=[
                        {"role": "system", "content": "You are a helpful assistant that validates bioinformatics analysis results."},
                        {"role": "user", "content": prompt}
                    ],
                    temperature=0.2
                )
                return response.choices[0].message.content
        except Exception as e:
            print(f"OpenAI API error: {e}")
            return None

    async def _call_local_model(self, prompt: str) -> Optional[str]:
        """调用本地模型
        
        Args:
            prompt: 提示词
            
        Returns:
            Optional[str]: 模型响应
        """
        try:
            import httpx
            
            async with asyncio.timeout(self.timeout):
                async with httpx.AsyncClient() as client:
                    response = await client.post(
                        f"{self.local_model_url}/generate",
                        json={"prompt": prompt, "max_tokens": 1000}
                    )
                    response.raise_for_status()
                    return response.json().get("text")
        except Exception as e:
            print(f"Local model error: {e}")
            return None

    def _parse_llm_response(self, response: str) -> Optional[ValidationResult]:
        """解析LLM响应为ValidationResult
        
        Args:
            response: LLM响应
            
        Returns:
            Optional[ValidationResult]: 解析后的验证结果
        """
        try:
            # 提取JSON部分
            import re
            json_match = re.search(r'\{[^}]*\}', response)
            if not json_match:
                return None
            
            json_str = json_match.group(0)
            data = json.loads(json_str)
            
            # 验证必要字段
            if 'passed' not in data or 'score' not in data:
                return None
            
            return ValidationResult(
                passed=data['passed'],
                score=data['score'],
                issues=data.get('issues', []),
                suggestion=data.get('suggestion')
            )
        except Exception as e:
            print(f"Failed to parse LLM response: {e}")
            return None

    def _fallback_validation(self, content: str) -> ValidationResult:
        """降级验证（当LLM失败时使用）
        
        Args:
            content: 要验证的内容
            
        Returns:
            ValidationResult: 降级验证结果
        """
        # 简单的规则验证
        issues = []
        
        # 检查内容是否为空
        if not content or content.strip() == '':
            issues.append("输出内容为空")
        
        # 检查是否包含错误信息
        error_keywords = ['error', 'failed', 'exception', 'traceback']
        for keyword in error_keywords:
            if keyword.lower() in content.lower():
                issues.append(f"可能包含错误信息: {keyword}")
        
        passed = len(issues) == 0
        score = 1.0 if passed else 0.5
        
        return ValidationResult(
            passed=passed,
            score=score,
            issues=issues,
            suggestion="使用规则验证，建议使用LLM进行更详细的分析"
        )

    async def validate(self, content: str, task_name: str) -> ValidationResult:
        """验证内容
        
        Args:
            content: 要验证的内容
            task_name: 任务名称
            
        Returns:
            ValidationResult: 验证结果
        """
        # 检查缓存
        cached_result = self._get_from_cache(content, task_name)
        if cached_result:
            return cached_result
        
        # 构造Prompt
        prompt = f"分析以下{task_name}输出，检查常见错误：\n{content}\n\n返回JSON：{{\"passed\": bool, \"score\": float, \"issues\": [str], \"suggestion\": str}}"
        
        # 调用LLM
        if self.model_type == 'openai' and self.openai_api_key:
            response = await self._call_openai(prompt)
        else:
            response = await self._call_local_model(prompt)
        
        # 解析响应
        result = self._parse_llm_response(response) if response else None
        
        # 降级处理
        if not result:
            result = self._fallback_validation(content)
        
        # 保存到缓存
        self._save_to_cache(content, task_name, result)
        
        return result


class RuleValidator:
    def __init__(self, config: dict = None):
        """初始化规则验证器
        
        Args:
            config: 配置字典，包含验证阈值等
        """
        self.config = config or {}
        self.min_file_size = self.config.get('min_file_size', 10)  # 最小文件大小（字节）
        self.min_peaks = self.config.get('min_peaks', 100)  # 最小峰值数量
        self.max_log2fc = self.config.get('max_log2fc', 10)  # 最大log2FC值
        self.min_log2fc_std = self.config.get('min_log2fc_std', 0.1)  # 最小log2FC标准差

    def check_file_not_empty(self, filepath: str) -> ValidationResult:
        """检查文件是否非空
        
        Args:
            filepath: 文件路径
            
        Returns:
            ValidationResult: 验证结果
        """
        issues = []
        
        try:
            if not os.path.exists(filepath):
                issues.append(f"文件不存在: {filepath}")
            elif os.path.getsize(filepath) < self.min_file_size:
                issues.append(f"文件为空或太小（小于{self.min_file_size}字节）")
        except Exception as e:
            issues.append(f"检查文件时出错: {str(e)}")
        
        passed = len(issues) == 0
        score = 1.0 if passed else 0.0
        
        return ValidationResult(
            passed=passed,
            score=score,
            issues=issues,
            suggestion="确保文件存在且包含有效内容"
        )

    def check_csv_columns(self, filepath: str, required_columns: List[str]) -> ValidationResult:
        """检查CSV文件列完整性
        
        Args:
            filepath: CSV文件路径
            required_columns: 必需的列名列表
            
        Returns:
            ValidationResult: 验证结果
        """
        issues = []
        
        try:
            # 首先检查文件是否存在
            file_check = self.check_file_not_empty(filepath)
            if not file_check.passed:
                return file_check
            
            # 读取CSV文件
            df = pd.read_csv(filepath)
            
            # 检查必需列
            missing_columns = [col for col in required_columns if col not in df.columns]
            if missing_columns:
                issues.append(f"缺少必需列: {', '.join(missing_columns)}")
            
            # 检查数据行数
            if len(df) == 0:
                issues.append("CSV文件为空（无数据行）")
        except Exception as e:
            issues.append(f"读取CSV文件时出错: {str(e)}")
        
        passed = len(issues) == 0
        score = 1.0 if passed else 0.5
        
        return ValidationResult(
            passed=passed,
            score=score,
            issues=issues,
            suggestion="确保CSV文件包含所有必需的列且有数据"
        )

    def check_log2fc_distribution(self, filepath: str) -> ValidationResult:
        """检查差异表达结果分布合理性
        
        Args:
            filepath: 差异表达结果文件路径
            
        Returns:
            ValidationResult: 验证结果
        """
        issues = []
        
        try:
            # 首先检查文件是否存在
            file_check = self.check_file_not_empty(filepath)
            if not file_check.passed:
                return file_check
            
            # 读取文件
            df = pd.read_csv(filepath)
            
            # 检查是否有log2FC列
            log2fc_columns = [col for col in df.columns if 'log2' in col.lower() or 'fc' in col.lower()]
            if not log2fc_columns:
                issues.append("未找到log2FC列")
            else:
                # 使用第一个log2FC列
                log2fc_col = log2fc_columns[0]
                log2fc_values = df[log2fc_col].dropna()
                
                # 检查值范围
                if len(log2fc_values) > 0:
                    max_abs_log2fc = np.max(np.abs(log2fc_values))
                    if max_abs_log2fc > self.max_log2fc:
                        issues.append(f"log2FC值过大（最大值: {max_abs_log2fc}）")
                    
                    # 检查分布
                    std_log2fc = np.std(log2fc_values)
                    if std_log2fc < self.min_log2fc_std:
                        issues.append(f"log2FC分布过于集中（标准差: {std_log2fc}）")
                    
                    # 检查是否有异常值
                    q1 = np.percentile(log2fc_values, 25)
                    q3 = np.percentile(log2fc_values, 75)
                    iqr = q3 - q1
                    outliers = log2fc_values[(log2fc_values < q1 - 1.5 * iqr) | (log2fc_values > q3 + 1.5 * iqr)]
                    if len(outliers) > len(log2fc_values) * 0.1:  # 超过10%的异常值
                        issues.append("存在过多的log2FC异常值")
        except Exception as e:
            issues.append(f"检查log2FC分布时出错: {str(e)}")
        
        passed = len(issues) == 0
        score = 1.0 if passed else 0.6
        
        return ValidationResult(
            passed=passed,
            score=score,
            issues=issues,
            suggestion="确保差异表达结果的log2FC分布合理，没有极端异常值"
        )

    def check_bam_index(self, bam_path: str) -> ValidationResult:
        """检查BAM文件索引是否存在
        
        Args:
            bam_path: BAM文件路径
            
        Returns:
            ValidationResult: 验证结果
        """
        issues = []
        
        try:
            # 首先检查BAM文件是否存在
            file_check = self.check_file_not_empty(bam_path)
            if not file_check.passed:
                return file_check
            
            # 检查索引文件
            index_path = bam_path + '.bai'
            if not os.path.exists(index_path):
                issues.append(f"BAM索引文件不存在: {index_path}")
            elif os.path.getsize(index_path) < self.min_file_size:
                issues.append("BAM索引文件为空或太小")
        except Exception as e:
            issues.append(f"检查BAM索引时出错: {str(e)}")
        
        passed = len(issues) == 0
        score = 1.0 if passed else 0.0
        
        return ValidationResult(
            passed=passed,
            score=score,
            issues=issues,
            suggestion="确保BAM文件有对应的索引文件（.bai）"
        )

    def check_peak_count(self, bed_path: str, min_peaks: int = None) -> ValidationResult:
        """检查峰值数量是否合理
        
        Args:
            bed_path: BED文件路径
            min_peaks: 最小峰值数量，默认为配置中的值
            
        Returns:
            ValidationResult: 验证结果
        """
        if min_peaks is None:
            min_peaks = self.min_peaks
        
        issues = []
        
        try:
            # 首先检查文件是否存在
            file_check = self.check_file_not_empty(bed_path)
            if not file_check.passed:
                return file_check
            
            # 计算峰值数量
            with open(bed_path, 'r') as f:
                peak_count = sum(1 for line in f if line.strip() and not line.startswith('#'))
            
            if peak_count < min_peaks:
                issues.append(f"峰值数量过少（{peak_count}个，最少需要{min_peaks}个）")
            elif peak_count > 1000000:  # 峰值数量过多
                issues.append(f"峰值数量过多（{peak_count}个），可能存在问题")
        except Exception as e:
            issues.append(f"检查峰值数量时出错: {str(e)}")
        
        passed = len(issues) == 0
        score = 1.0 if passed else 0.5
        
        return ValidationResult(
            passed=passed,
            score=score,
            issues=issues,
            suggestion="确保峰值数量在合理范围内"
        )
