import pytest
import os
import tempfile
from unittest import mock
import asyncio
from ..src.models import UserRequest, TaskType, PipelineTemplate, PipelineStep
from ..src.planner import Planner
from ..src.executor import AsyncPipelineExecutor
from ..src.arbiter import Arbiter, AdaptiveIterator
from ..src.registry import PipelineRegistry
from ..src.webui.app import AutoBAWebUI


def test_rnaseq_pipeline_mock(temp_dir):
    """使用mock数据测试完整RNA-seq流程"""
    # 创建mock数据文件
    mock_expression_path = os.path.join(temp_dir, "mock_expression.csv")
    with open(mock_expression_path, "w") as f:
        f.write("gene,log2FoldChange,pvalue\n")
        f.write("TP53,2.5,0.001\n")
        f.write("BRCA1,1.8,0.005\n")
        f.write("EGFR,-1.2,0.01\n")
    
    # 创建测试配置
    config = {
        "llm": {
            "mode": "mock"
        },
        "executor": {
            "timeout": 3600,
            "command_whitelist": []
        },
        "output_dir": temp_dir
    }
    
    # 创建执行器
    executor = AsyncPipelineExecutor(config.get('executor', {}))
    
    # 创建测试步骤
    steps = [
        PipelineStep(
            step_id=1,
            tool_name="echo",
            tool_version="1.0",
            command_template="echo 'Processing {input_file}' > {output_file}",
            input_files=["{input_file}"],
            output_files=["{output_file}"],
            parameters={}
        )
    ]
    
    # 创建测试模板
    template = PipelineTemplate(
        name="test_rnaseq",
        description="Test RNA-seq pipeline",
        task_type="deterministic",
        steps=steps,
        required_inputs=["input_file"],
        output_patterns=["{output_file}"]
    )
    
    # 创建用户请求
    user_request = UserRequest(
        data_paths=[mock_expression_path],
        data_description="RNA-seq expression data",
        goal="Differential expression analysis",
        task_type=TaskType.DETERMINISTIC
    )
    
    # 执行流程
    async def run_pipeline():
        context = {
            "input_file": mock_expression_path,
            "output_file": os.path.join(temp_dir, "output.txt")
        }
        records = await executor.run_pipeline(
            template=template,
            user_input=user_request,
            workdir=temp_dir
        )
        return records
    
    records = asyncio.run(run_pipeline())
    assert len(records) == 1
    assert records[0].status == "success"

@pytest.mark.asyncio
async def test_llm_planner_integration(temp_dir):
    """测试规划器与执行器集成"""
    # 创建测试配置
    config = {
        "llm": {
            "mode": "mock"
        },
        "executor": {
            "timeout": 3600,
            "command_whitelist": []
        },
        "output_dir": temp_dir
    }
    
    # 创建规划器和执行器
    planner = Planner(config)
    executor = AsyncPipelineExecutor(config.get('executor', {}))
    
    # 创建用户请求
    user_request = UserRequest(
        data_paths=["sample1.fastq", "sample2.fastq"],
        data_description="RNA-seq data",
        goal="Differential expression analysis",
        task_type=TaskType.DETERMINISTIC
    )
    
    # 生成计划
    template = await planner.plan(user_request)
    assert template is not None
    assert template.name is not None
    
    # 执行计划（使用mock）
    with mock.patch('asyncio.create_subprocess_shell') as mock_create_subprocess:
        # 配置mock
        mock_process = mock.Mock()
        mock_process.communicate.return_value = (b"Command executed successfully", b"")
        mock_process.returncode = 0
        mock_create_subprocess.return_value = mock_process
        
        # 执行流程
        context = {
            "output_dir": temp_dir
        }
        records = await executor.run_pipeline(
            template=template,
            user_input=user_request,
            workdir=temp_dir
        )
        
        assert len(records) > 0
        assert records[0].status == "success"


def test_arbiter_retry_logic():
    """测试失败重试机制"""
    # 创建测试配置
    config = {
        "max_retries": 3,
        "score_improvement_threshold": 0.1
    }
    
    # 创建仲裁器
    arbiter = Arbiter(config)
    
    # 模拟执行历史
    history = []
    
    # 第一次验证：失败，分数低
    validation1 = mock.Mock(
        passed=False,
        score=0.5
    )
    decision1 = arbiter.decide_next_action(history, validation1)
    assert decision1 == "retry"
    history.append({"status": "failed", "score": 0.5})
    
    # 第二次验证：失败，分数有所提高
    validation2 = mock.Mock(
        passed=False,
        score=0.65
    )
    decision2 = arbiter.decide_next_action(history, validation2)
    assert decision2 == "retry"
    history.append({"status": "failed", "score": 0.65})
    
    # 第三次验证：失败，分数提高不足
    validation3 = mock.Mock(
        passed=False,
        score=0.7
    )
    decision3 = arbiter.decide_next_action(history, validation3)
    assert decision3 == "fallback"


def test_webui_api():
    """测试Gradio接口"""
    # 创建测试配置
    config = {
        "llm": {
            "mode": "mock"
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
    
    # 创建WebUI实例
    webui = AutoBAWebUI(config)
    assert webui is not None
    
    # 测试get_pipeline_list方法
    pipeline_list = webui.get_pipeline_list()
    assert isinstance(pipeline_list, list)
    
    # 测试get_task_history方法
    task_history = webui.get_task_history()
    assert isinstance(task_history, str)
    assert "任务历史" in task_history


def test_adaptive_iterator_integration():
    """测试AdaptiveIterator集成"""
    # 创建测试配置
    config = {
        "param_space": {
            "learning_rate": {
                "type": "continuous",
                "min": 0.001,
                "max": 0.1,
                "default": 0.01
            },
            "batch_size": {
                "type": "discrete",
                "values": [16, 32, 64],
                "default": 32
            }
        }
    }
    
    # 创建自适应迭代器
    iterator = AdaptiveIterator(config)
    
    # 模拟执行记录
    for i in range(3):
        # 模拟不同的验证结果
        score = 0.5 + i * 0.1
        record = {
            "params": {
                "learning_rate": 0.01,
                "batch_size": 32
            },
            "validation": {
                "score": score,
                "passed": score > 0.7
            }
        }
        iterator.add_execution_record(record)
    
    # 获取下一组参数
    next_params = iterator.get_next_params()
    assert "learning_rate" in next_params
    assert "batch_size" in next_params
    assert 0.001 <= next_params["learning_rate"] <= 0.1
    assert next_params["batch_size"] in [16, 32, 64]
