import pytest
import os
import tempfile
import yaml
from ..src.registry import PipelineRegistry
from ..src.models import PipelineTemplate


def test_pipeline_registry_initialization():
    """测试PipelineRegistry初始化"""
    registry = PipelineRegistry()
    assert registry is not None


def test_pipeline_registry_discover(temp_dir):
    """测试PipelineRegistry discover方法"""
    # 创建测试配置文件
    pipelines_dir = os.path.join(temp_dir, "pipelines")
    os.makedirs(pipelines_dir, exist_ok=True)
    
    # 创建测试YAML文件
    test_yaml = {
        "name": "test_pipeline",
        "description": "Test pipeline",
        "task_type": "deterministic",
        "required_inputs": ["input1", "input2"],
        "output_patterns": ["output1.txt"],
        "steps": [
            {
                "step_id": 1,
                "tool_name": "echo",
                "tool_version": "1.0",
                "command_template": "echo {input1} > {output1}",
                "input_files": ["{input1}"],
                "output_files": ["{output1}"],
                "parameters": {}
            }
        ]
    }
    
    yaml_path = os.path.join(pipelines_dir, "test_pipeline.yaml")
    with open(yaml_path, "w") as f:
        yaml.dump(test_yaml, f)
    
    # 测试discover方法
    registry = PipelineRegistry(pipelines_dir=pipelines_dir)
    templates = registry.discover()
    assert len(templates) == 1
    assert templates[0].name == "test_pipeline"


def test_pipeline_registry_get(temp_dir):
    """测试PipelineRegistry get方法"""
    # 创建测试配置文件
    pipelines_dir = os.path.join(temp_dir, "pipelines")
    os.makedirs(pipelines_dir, exist_ok=True)
    
    # 创建测试YAML文件
    test_yaml = {
        "name": "test_pipeline",
        "description": "Test pipeline",
        "task_type": "deterministic",
        "required_inputs": ["input1", "input2"],
        "output_patterns": ["output1.txt"],
        "steps": [
            {
                "step_id": 1,
                "tool_name": "echo",
                "tool_version": "1.0",
                "command_template": "echo {input1} > {output1}",
                "input_files": ["{input1}"],
                "output_files": ["{output1}"],
                "parameters": {}
            }
        ]
    }
    
    yaml_path = os.path.join(pipelines_dir, "test_pipeline.yaml")
    with open(yaml_path, "w") as f:
        yaml.dump(test_yaml, f)
    
    # 测试get方法
    registry = PipelineRegistry(pipelines_dir=pipelines_dir)
    template = registry.get("test_pipeline")
    assert template is not None
    assert template.name == "test_pipeline"
    
    # 测试获取不存在的模板
    template_none = registry.get("non_existent")
    assert template_none is None


def test_pipeline_registry_validate_template():
    """测试PipelineRegistry validate_template方法"""
    # 创建一个有效的模板
    from ..src.models import PipelineStep
    steps = [
        PipelineStep(
            step_id=1,
            tool_name="echo",
            tool_version="1.0",
            command_template="echo {input1} > {output1}",
            input_files=["{input1}"],
            output_files=["{output1}"],
            parameters={}
        )
    ]
    valid_template = PipelineTemplate(
        name="valid_pipeline",
        description="Valid pipeline",
        task_type="deterministic",
        steps=steps,
        required_inputs=["input1"],
        output_patterns=["output1.txt"]
    )
    
    registry = PipelineRegistry()
    assert registry.validate_template(valid_template) is True


def test_pipeline_registry_get_dependencies():
    """测试PipelineRegistry get_dependencies方法"""
    # 创建一个包含工具的模板
    from ..src.models import PipelineStep
    steps = [
        PipelineStep(
            step_id=1,
            tool_name="echo",
            tool_version="1.0",
            command_template="echo {input1} > {output1}",
            input_files=["{input1}"],
            output_files=["{output1}"],
            parameters={}
        ),
        PipelineStep(
            step_id=2,
            tool_name="cat",
            tool_version="1.0",
            command_template="cat {output1} > {output2}",
            input_files=["{output1}"],
            output_files=["{output2}"],
            parameters={}
        )
    ]
    template = PipelineTemplate(
        name="test_pipeline",
        description="Test pipeline",
        task_type="deterministic",
        steps=steps,
        required_inputs=["input1"],
        output_patterns=["output2.txt"]
    )
    
    registry = PipelineRegistry()
    dependencies = registry.get_dependencies(template)
    assert "echo" in dependencies
    assert "cat" in dependencies
    assert len(dependencies) == 2


def test_pipeline_registry_to_json_schema():
    """测试PipelineRegistry to_json_schema方法"""
    # 创建一个模板
    from ..src.models import PipelineStep
    steps = [
        PipelineStep(
            step_id=1,
            tool_name="echo",
            tool_version="1.0",
            command_template="echo {input1} > {output1}",
            input_files=["{input1}"],
            output_files=["{output1}"],
            parameters={
                "param1": "default_value"
            }
        )
    ]
    template = PipelineTemplate(
        name="test_pipeline",
        description="Test pipeline",
        task_type="deterministic",
        steps=steps,
        required_inputs=["input1"],
        output_patterns=["output1.txt"]
    )
    
    registry = PipelineRegistry()
    schema = registry.to_json_schema(template)
    assert schema is not None
    assert "type" in schema
    assert schema["type"] == "object"
    assert "properties" in schema


def test_pipeline_registry_template_inheritance(temp_dir):
    """测试PipelineRegistry模板继承"""
    # 创建测试配置文件
    pipelines_dir = os.path.join(temp_dir, "pipelines")
    os.makedirs(pipelines_dir, exist_ok=True)
    
    # 创建基础模板
    base_yaml = {
        "name": "base_pipeline",
        "description": "Base pipeline",
        "task_type": "deterministic",
        "required_inputs": ["input1"],
        "output_patterns": ["output1.txt"],
        "steps": [
            {
                "step_id": 1,
                "tool_name": "echo",
                "tool_version": "1.0",
                "command_template": "echo {input1} > {output1}",
                "input_files": ["{input1}"],
                "output_files": ["{output1}"],
                "parameters": {}
            }
        ]
    }
    
    # 创建继承模板
    child_yaml = {
        "name": "child_pipeline",
        "extends": "base_pipeline",
        "description": "Child pipeline",
        "required_inputs": ["input1", "input2"],
        "steps": [
            {
                "step_id": 2,
                "tool_name": "cat",
                "tool_version": "1.0",
                "command_template": "cat {output1} {input2} > {output2}",
                "input_files": ["{output1}", "{input2}"],
                "output_files": ["{output2}"],
                "parameters": {}
            }
        ]
    }
    
    # 写入文件
    base_path = os.path.join(pipelines_dir, "base_pipeline.yaml")
    with open(base_path, "w") as f:
        yaml.dump(base_yaml, f)
    
    child_path = os.path.join(pipelines_dir, "child_pipeline.yaml")
    with open(child_path, "w") as f:
        yaml.dump(child_yaml, f)
    
    # 测试模板继承
    registry = PipelineRegistry(pipelines_dir=pipelines_dir)
    child_template = registry.get("child_pipeline")
    assert child_template is not None
    assert child_template.name == "child_pipeline"
    assert len(child_template.steps) == 2  # 应该包含基础模板的步骤
    assert len(child_template.required_inputs) == 2  # 应该包含新增的输入
