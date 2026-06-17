import pytest
import os
import tempfile
from ..src.models import PipelineStep, PipelineTemplate, UserRequest, TaskType


@pytest.fixture
def temp_dir():
    """临时目录fixture"""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def sample_pipeline_step():
    """示例PipelineStep"""
    return PipelineStep(
        step_id=1,
        tool_name="Trimmomatic",
        tool_version="0.39",
        command_template="trimmomatic PE {input1} {input2} {output1} {output2} ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
        input_files=["{input1}", "{input2}"],
        output_files=["{output1}", "{output2}"],
        parameters={
            "adapter": "adapters.fa",
            "leading": 3,
            "trailing": 3,
            "slidingwindow": "4:15",
            "minlen": 36
        }
    )


@pytest.fixture
def sample_pipeline_template():
    """示例PipelineTemplate"""
    steps = [
        PipelineStep(
            step_id=1,
            tool_name="Trimmomatic",
            tool_version="0.39",
            command_template="trimmomatic PE {input1} {input2} {output1} {output2} ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
            input_files=["{input1}", "{input2}"],
            output_files=["{output1}", "{output2}"],
            parameters={
                "adapter": "adapters.fa",
                "leading": 3,
                "trailing": 3,
                "slidingwindow": "4:15",
                "minlen": 36
            }
        ),
        PipelineStep(
            step_id=2,
            tool_name="HISAT2",
            tool_version="2.2.1",
            command_template="hisat2 -x {genome_index} -1 {input1} -2 {input2} -S {output}",
            input_files=["{input1}", "{input2}"],
            output_files=["{output}"],
            parameters={
                "genome_index": "genome_index"
            }
        )
    ]
    return PipelineTemplate(
        name="rnaseq_diffexpr",
        description="RNA-seq差异表达分析流程",
        task_type="deterministic",
        steps=steps,
        required_inputs=["fastq_files", "genome_index", "annotation_gtf"],
        output_patterns=["deseq2_results.csv"]
    )


@pytest.fixture
def sample_user_request():
    """示例UserRequest"""
    return UserRequest(
        data_paths=["sample1.fastq", "sample2.fastq"],
        data_description="RNA-seq数据，两个样本",
        goal="差异表达分析",
        task_type=TaskType.DETERMINISTIC
    )


@pytest.fixture
def sample_config():
    """示例配置"""
    return {
        "llm": {
            "mode": "mock",
            "mock": {}
        },
        "executor": {
            "timeout": 3600,
            "command_whitelist": []
        },
        "xai": {
            "language": "zh"
        },
        "output_dir": "./output"
    }
