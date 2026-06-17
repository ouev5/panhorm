from typing import List, Dict, Any, Optional
from models import ValidationResult, ExecutionRecord


class Arbiter:
    def __init__(self, config: dict = None):
        """初始化仲裁器
        
        Args:
            config: 配置字典，包含决策阈值等
        """
        self.config = config or {}
        self.max_retries = self.config.get('max_retries', 3)  # 最大重试次数
        self.score_improvement_threshold = self.config.get('score_improvement_threshold', 0.1)  # 分数提升阈值

    def decide_next_action(self, history: List[ExecutionRecord], last_validation: ValidationResult) -> str:
        """决定下一步行动
        
        Args:
            history: 执行历史记录
            last_validation: 最后一次验证结果
            
        Returns:
            str: 决策结果，可选值: "retry", "fallback", "success", "abort"
        """
        # 规则1: 验证通过 → success
        if last_validation.passed:
            return "success"
        
        # 规则2: 重试次数达到最大值 → abort
        retry_count = len([record for record in history if record.status == "RETRY"])
        if retry_count >= self.max_retries:
            return "abort"
        
        # 规则3: 验证分数相比上次提升小于阈值 → fallback
        if len(history) >= 2:
            # 找到上一次的验证结果
            previous_validation = None
            for record in reversed(history[:-1]):
                if hasattr(record, 'validation_result') and record.validation_result:
                    previous_validation = record.validation_result
                    break
            
            if previous_validation:
                score_improvement = last_validation.score - previous_validation.score
                if score_improvement < self.score_improvement_threshold:
                    return "fallback"
        
        # 规则4: 否则 → retry
        return "retry"


class AdaptiveIterator:
    def __init__(self, config: dict = None):
        """初始化自适应迭代器
        
        Args:
            config: 配置字典，包含参数搜索空间等
        """
        self.config = config or {}
        self.execution_history = []  # 执行历史
        self.param_search_space = self.config.get('param_search_space', {})  # 参数搜索空间
        self.current_params = self.config.get('initial_params', {})  # 当前参数
        self.arbiter = Arbiter(self.config.get('arbiter_config', {}))  # 仲裁器

    def add_execution_record(self, record: ExecutionRecord):
        """添加执行记录
        
        Args:
            record: 执行记录
        """
        self.execution_history.append(record)

    def get_next_params(self, last_validation: ValidationResult) -> Dict[str, Any]:
        """获取下一组参数
        
        Args:
            last_validation: 最后一次验证结果
            
        Returns:
            Dict[str, Any]: 下一组参数
        """
        # 决定下一步行动
        action = self.arbiter.decide_next_action(self.execution_history, last_validation)
        
        if action == "success":
            # 成功，返回当前参数
            return self.current_params
        elif action == "abort":
            # 终止，返回None
            return None
        elif action == "fallback":
            # 回退，使用默认参数
            return self.config.get('default_params', {})
        elif action == "retry":
            # 重试，调整参数
            return self._adjust_params(last_validation)
        
        return self.current_params

    def _adjust_params(self, validation: ValidationResult) -> Dict[str, Any]:
        """调整参数
        
        Args:
            validation: 验证结果
            
        Returns:
            Dict[str, Any]: 调整后的参数
        """
        new_params = self.current_params.copy()
        
        # 根据验证结果和参数搜索空间调整参数
        for param_name, param_space in self.param_search_space.items():
            if param_name in new_params:
                current_value = new_params[param_name]
                
                # 检查参数类型并调整
                if isinstance(param_space, dict):
                    # 连续参数空间
                    if 'min' in param_space and 'max' in param_space and 'step' in param_space:
                        # 根据验证分数调整
                        if validation.score < 0.5:
                            # 分数低，向默认值靠拢
                            default_value = param_space.get('default', (param_space['min'] + param_space['max']) / 2)
                            new_params[param_name] = self._lerp(current_value, default_value, 0.5)
                        else:
                            # 分数高，继续优化
                            # 这里可以实现更复杂的优化策略
                            pass
                elif isinstance(param_space, list):
                    # 离散参数空间
                    if current_value in param_space:
                        current_idx = param_space.index(current_value)
                        # 简单的轮询策略
                        next_idx = (current_idx + 1) % len(param_space)
                        new_params[param_name] = param_space[next_idx]
        
        self.current_params = new_params
        return new_params

    def _lerp(self, a: float, b: float, t: float) -> float:
        """线性插值
        
        Args:
            a: 起始值
            b: 结束值
            t: 插值参数 (0-1)
            
        Returns:
            float: 插值结果
        """
        return a + (b - a) * t

    def define_param_space(self, param_name: str, space: Any):
        """定义参数搜索空间
        
        Args:
            param_name: 参数名称
            space: 参数搜索空间，可以是列表（离散值）或字典（连续值，包含min, max, step, default）
        """
        self.param_search_space[param_name] = space

    def get_execution_history(self) -> List[ExecutionRecord]:
        """获取执行历史
        
        Returns:
            List[ExecutionRecord]: 执行历史记录
        """
        return self.execution_history
