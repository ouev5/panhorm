from pydantic import BaseModel, Field, field_validator, model_validator
from typing import List, Dict, Any, Optional
from enum import Enum
from datetime import datetime


class TaskType(str, Enum):
    DETERMINISTIC = "deterministic"
    AMBIGUOUS = "ambiguous"


class ExecutionStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    SUCCESS = "success"
    FAILED = "failed"
    RETRY = "retry"
    ABORTED = "aborted"


class PipelineStep(BaseModel):
    step_id: int = Field(..., ge=1, description="步骤ID，必须大于0")
    tool_name: str = Field(..., description="工具名称")
    tool_version: str = Field(..., description="工具版本")
    command_template: str = Field(..., description="命令模板")
    input_files: List[str] = Field(default_factory=list, description="输入文件列表")
    output_files: List[str] = Field(default_factory=list, description="输出文件列表")
    parameters: Dict[str, Any] = Field(default_factory=dict, description="参数字典")

    class Config:
        frozen = True

    @field_validator('command_template')
    def validate_command_template(cls, v):
        # 检查危险操作
        dangerous_patterns = [
            'rm -rf',
            'curl.*\\|.*bash',
            'wget.*\\|.*bash',
            'chmod.*777',
            'sudo ',
            'su -'
        ]
        for pattern in dangerous_patterns:
            if pattern in v:
                raise ValueError(f"命令模板包含危险操作: {pattern}")
        return v

    @field_validator('step_id')
    def validate_step_id(cls, v):
        if not 1 <= v <= 1000:
            raise ValueError("step_id必须在1-1000范围内")
        return v

    @model_validator(mode='after')
    def validate_files(self):
        if not self.input_files:
            raise ValueError("input_files不能为空")
        if not self.output_files:
            raise ValueError("output_files不能为空")
        return self

    def render_command(self, context: dict) -> str:
        """安全地替换命令模板中的变量"""
        command = self.command_template
        # 安全地替换变量，防止注入
        for key, value in context.items():
            # 确保值是字符串
            safe_value = str(value)
            # 避免危险字符
            safe_value = safe_value.replace(';', '').replace('|', '').replace('&', '')
            # 替换模板变量
            command = command.replace(f"{{{{{key}}}}}", safe_value)
        return command

    def to_dict(self) -> Dict[str, Any]:
        return self.model_dump()

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "PipelineStep":
        return cls(**data)


class PipelineTemplate(BaseModel):
    name: str = Field(..., description="管道模板名称")
    description: str = Field(..., description="管道模板描述")
    task_type: TaskType = Field(..., description="任务类型")
    steps: List[PipelineStep] = Field(default_factory=list, description="管道步骤列表")
    required_inputs: List[str] = Field(default_factory=list, description="必需的输入文件列表")
    output_patterns: List[str] = Field(default_factory=list, description="输出文件模式列表")

    class Config:
        frozen = True

    def to_dict(self) -> Dict[str, Any]:
        return self.model_dump()

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "PipelineTemplate":
        return cls(**data)


class UserRequest(BaseModel):
    data_paths: List[str] = Field(..., description="数据文件路径列表")
    data_description: str = Field(..., description="数据描述")
    goal: str = Field(..., description="分析目标")
    task_type: TaskType = Field(..., description="任务类型")

    def to_dict(self) -> Dict[str, Any]:
        return self.model_dump()

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "UserRequest":
        return cls(**data)


class ExecutionRecord(BaseModel):
    execution_id: str = Field(..., description="执行ID")
    pipeline_name: str = Field(..., description="管道名称")
    start_time: datetime = Field(default_factory=datetime.now, description="开始时间")
    end_time: Optional[datetime] = Field(None, description="结束时间")
    status: ExecutionStatus = Field(default=ExecutionStatus.PENDING, description="执行状态")
    step_executions: List[Dict[str, Any]] = Field(default_factory=list, description="步骤执行记录")
    error_message: Optional[str] = Field(None, description="错误信息")

    def to_dict(self) -> Dict[str, Any]:
        return self.model_dump()

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ExecutionRecord":
        return cls(**data)


class ValidationResult(BaseModel):
    passed: bool = Field(..., description="验证是否通过")
    score: float = Field(..., ge=0.0, le=1.0, description="验证分数，范围0-1")
    issues: List[str] = Field(default_factory=list, description="问题列表")
    suggestion: Optional[str] = Field(None, description="改进建议")

    class Config:
        frozen = True

    def to_dict(self) -> Dict[str, Any]:
        return self.model_dump()

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ValidationResult":
        return cls(**data)


class XAIOutput(BaseModel):
    confidence: float = Field(..., ge=0.0, le=1.0, description="置信度，范围0-1")
    evidence_sources: List[str] = Field(default_factory=list, description="证据来源列表")
    feature_contributions: Dict[str, float] = Field(default_factory=dict, description="特征贡献度")
    natural_language_explanation: str = Field(..., description="自然语言解释")

    class Config:
        frozen = True

    def to_dict(self) -> Dict[str, Any]:
        return self.model_dump()

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "XAIOutput":
        return cls(**data)


class FinalResult(BaseModel):
    status: ExecutionStatus = Field(..., description="执行状态")
    execution_summary: str = Field(..., description="执行摘要")
    output_files: List[str] = Field(default_factory=list, description="输出文件列表")
    xai_report: Optional[XAIOutput] = Field(None, description="XAI分析报告")
    logs: List[str] = Field(default_factory=list, description="日志列表")

    class Config:
        frozen = True

    def to_dict(self) -> Dict[str, Any]:
        return self.model_dump()

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "FinalResult":
        return cls(**data)
