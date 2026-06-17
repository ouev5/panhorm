import os
import yaml
from typing import Dict, List, Optional, Any
from models import PipelineTemplate, PipelineStep, TaskType


class PipelineRegistry:
    def __init__(self, pipelines_dir: str = "config/pipelines"):
        """初始化管道注册表
        
        Args:
            pipelines_dir: 管道模板目录
        """
        self.pipelines_dir = pipelines_dir
        self.templates: Dict[str, PipelineTemplate] = {}
        self._load_templates()

    def _load_templates(self):
        """加载所有模板"""
        if not os.path.exists(self.pipelines_dir):
            return
        
        # 先加载所有基础模板
        yaml_files = [f for f in os.listdir(self.pipelines_dir) if f.endswith('.yaml') or f.endswith('.yml')]
        
        # 分两步加载：先加载所有模板，再处理继承
        raw_templates = {}
        for yaml_file in yaml_files:
            file_path = os.path.join(self.pipelines_dir, yaml_file)
            with open(file_path, 'r', encoding='utf-8') as f:
                try:
                    raw_template = yaml.safe_load(f)
                    if raw_template and 'name' in raw_template:
                        raw_templates[raw_template['name']] = raw_template
                except yaml.YAMLError as e:
                    print(f"Error loading {yaml_file}: {e}")
        
        # 处理模板继承并转换为PipelineTemplate对象
        for name, raw_template in raw_templates.items():
            self._process_template(name, raw_template, raw_templates)

    def _process_template(self, name: str, raw_template: Dict[str, Any], raw_templates: Dict[str, Dict[str, Any]]):
        """处理模板继承并转换为PipelineTemplate对象
        
        Args:
            name: 模板名称
            raw_template: 原始模板数据
            raw_templates: 所有原始模板
        """
        # 处理继承
        if 'extends' in raw_template:
            base_name = raw_template['extends']
            if base_name in raw_templates:
                # 复制基础模板
                base_template = raw_templates[base_name].copy()
                # 合并当前模板的字段
                for key, value in raw_template.items():
                    if key == 'steps' and 'steps' in base_template:
                        # 合并步骤
                        base_template['steps'].extend(value)
                    else:
                        base_template[key] = value
                raw_template = base_template
        
        # 转换为PipelineTemplate对象
        try:
            # 转换步骤
            steps = []
            if 'steps' in raw_template:
                for step_data in raw_template['steps']:
                    step = PipelineStep(**step_data)
                    steps.append(step)
            
            template = PipelineTemplate(
                name=raw_template['name'],
                description=raw_template.get('description', ''),
                task_type=TaskType(raw_template.get('task_type', 'deterministic')),
                steps=steps,
                required_inputs=raw_template.get('required_inputs', []),
                output_patterns=raw_template.get('output_patterns', [])
            )
            
            self.templates[name] = template
        except Exception as e:
            print(f"Error processing template {name}: {e}")

    def discover(self) -> List[PipelineTemplate]:
        """扫描并加载所有模板
        
        Returns:
            List[PipelineTemplate]: 所有加载的模板
        """
        self._load_templates()
        return list(self.templates.values())

    def get(self, name: str) -> Optional[PipelineTemplate]:
        """根据名称获取模板
        
        Args:
            name: 模板名称
            
        Returns:
            Optional[PipelineTemplate]: 模板对象，如果不存在则返回None
        """
        return self.templates.get(name)

    def validate_template(self, template: PipelineTemplate) -> bool:
        """验证模板完整性
        
        Args:
            template: 模板对象
            
        Returns:
            bool: 模板是否有效
        """
        # 检查必需字段
        if not template.name:
            return False
        
        # 检查步骤顺序
        step_ids = [step.step_id for step in template.steps]
        if sorted(step_ids) != step_ids:
            return False
        
        # 检查步骤ID是否唯一
        if len(step_ids) != len(set(step_ids)):
            return False
        
        return True

    def get_dependencies(self, template: PipelineTemplate) -> List[str]:
        """提取所需的外部工具列表
        
        Args:
            template: 模板对象
            
        Returns:
            List[str]: 工具列表
        """
        dependencies = []
        for step in template.steps:
            tool_name = step.tool_name
            if tool_name not in dependencies:
                dependencies.append(tool_name)
        return dependencies

    def to_json_schema(self, template: PipelineTemplate) -> Dict[str, Any]:
        """生成前端动态表单的JSON Schema
        
        Args:
            template: 模板对象
            
        Returns:
            Dict[str, Any]: JSON Schema
        """
        schema = {
            "type": "object",
            "properties": {},
            "required": []
        }
        
        # 为每个必需输入添加字段
        for input_name in template.required_inputs:
            schema["properties"][input_name] = {
                "type": "string",
                "title": input_name,
                "description": f"Path to {input_name}"
            }
            schema["required"].append(input_name)
        
        # 为每个步骤的参数添加字段
        for step in template.steps:
            for param_name, param_value in step.parameters.items():
                param_key = f"{step.tool_name}_{param_name}"
                
                # 确定参数类型
                param_type = "string"
                if isinstance(param_value, bool):
                    param_type = "boolean"
                elif isinstance(param_value, (int, float)):
                    param_type = "number"
                
                schema["properties"][param_key] = {
                    "type": param_type,
                    "title": f"{step.tool_name} - {param_name}",
                    "default": param_value
                }
        
        return schema
