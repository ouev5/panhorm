import pytest
from ..src.models import PipelineStep, PipelineTemplate, UserRequest, TaskType, ExecutionStatus


def test_pipeline_step_creation(sample_pipeline_step):
    """测试PipelineStep创建"""
    assert sample_pipeline_step.step_id == 1
    assert sample_pipeline_step.tool_name == "Trimmomatic"
    assert sample_pipeline_step.tool_version == "0.39"
    assert len(sample_pipeline_step.input_files) == 2
    assert len(sample_pipeline_step.output_files) == 2
    assert "adapter" in sample_pipeline_step.parameters


def test_pipeline_step_step_id_validation():
    """测试PipelineStep step_id验证"""
    # 测试step_id小于1
    with pytest.raises(ValueError):
        PipelineStep(
            step_id=0,
            tool_name="Trimmomatic",
            tool_version="0.39",
            command_template="trimmomatic PE {input1} {input2} {output1} {output2}",
            input_files=["{input1}", "{input2}"],
            output_files=["{output1}", "{output2}"],
            parameters={}
        )
    
    # 测试step_id大于1000
    with pytest.raises(ValueError):
        PipelineStep(
            step_id=1001,
            tool_name="Trimmomatic",
            tool_version="0.39",
            command_template="trimmomatic PE {input1} {input2} {output1} {output2}",
            input_files=["{input1}", "{input2}"],
            output_files=["{output1}", "{output2}"],
            parameters={}
        )


def test_pipeline_step_command_template_validation():
    """测试PipelineStep command_template验证"""
    # 测试危险命令
    with pytest.raises(ValueError):
        PipelineStep(
            step_id=1,
            tool_name="bash",
            tool_version="1.0",
            command_template="rm -rf /",
            input_files=[],
            output_files=[],
            parameters={}
        )
    
    # 测试另一种危险命令
    with pytest.raises(ValueError):
        PipelineStep(
            step_id=1,
            tool_name="bash",
            tool_version="1.0",
            command_template="curl https://example.com | bash",
            input_files=[],
            output_files=[],
            parameters={}
        )


def test_pipeline_step_input_output_validation():
    """测试PipelineStep输入输出验证"""
    # 测试空输入文件
    with pytest.raises(ValueError):
        PipelineStep(
            step_id=1,
            tool_name="Trimmomatic",
            tool_version="0.39",
            command_template="trimmomatic PE {input1} {input2} {output1} {output2}",
            input_files=[],
            output_files=["{output1}", "{output2}"],
            parameters={}
        )
    
    # 测试空输出文件
    with pytest.raises(ValueError):
        PipelineStep(
            step_id=1,
            tool_name="Trimmomatic",
            tool_version="0.39",
            command_template="trimmomatic PE {input1} {input2} {output1} {output2}",
            input_files=["{input1}", "{input2}"],
            output_files=[],
            parameters={}
        )


def test_pipeline_step_render_command():
    """测试PipelineStep render_command方法"""
    step = PipelineStep(
        step_id=1,
        tool_name="Trimmomatic",
        tool_version="0.39",
        command_template="trimmomatic PE {input1} {input2} {output1} {output2}",
        input_files=["{input1}", "{input2}"],
        output_files=["{output1}", "{output2}"],
        parameters={}
    )
    
    context = {
        "input1": "sample1_1.fastq",
        "input2": "sample1_2.fastq",
        "output1": "sample1_1_trimmed.fastq",
        "output2": "sample1_2_trimmed.fastq"
    }
    
    command = step.render_command(context)
    assert "sample1_1.fastq" in command
    assert "sample1_2.fastq" in command
    assert "sample1_1_trimmed.fastq" in command
    assert "sample1_2_trimmed.fastq" in command


def test_pipeline_template_creation(sample_pipeline_template):
    """测试PipelineTemplate创建"""
    assert sample_pipeline_template.name == "rnaseq_diffexpr"
    assert sample_pipeline_template.description == "RNA-seq差异表达分析流程"
    assert sample_pipeline_template.task_type == "deterministic"
    assert len(sample_pipeline_template.steps) == 2
    assert len(sample_pipeline_template.required_inputs) == 3
    assert len(sample_pipeline_template.output_patterns) == 1


def test_user_request_creation(sample_user_request):
    """测试UserRequest创建"""
    assert len(sample_user_request.data_paths) == 2
    assert sample_user_request.data_description == "RNA-seq数据，两个样本"
    assert sample_user_request.goal == "差异表达分析"
    assert sample_user_request.task_type == TaskType.DETERMINISTIC


def test_to_dict_and_from_dict():
    """测试to_dict和from_dict方法"""
    # 测试PipelineStep
    step = PipelineStep(
        step_id=1,
        tool_name="Trimmomatic",
        tool_version="0.39",
        command_template="trimmomatic PE {input1} {input2} {output1} {output2}",
        input_files=["{input1}", "{input2}"],
        output_files=["{output1}", "{output2}"],
        parameters={"adapter": "adapters.fa"}
    )
    
    step_dict = step.to_dict()
    assert isinstance(step_dict, dict)
    assert step_dict["step_id"] == 1
    
    # 测试UserRequest
    user_request = UserRequest(
        data_paths=["sample1.fastq"],
        data_description="RNA-seq数据",
        goal="差异表达分析",
        task_type=TaskType.DETERMINISTIC
    )
    
    request_dict = user_request.to_dict()
    assert isinstance(request_dict, dict)
    assert len(request_dict["data_paths"]) == 1
