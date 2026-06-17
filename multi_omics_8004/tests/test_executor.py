import pytest
import os
import tempfile
from unittest import mock
from ..src.executor import PipelineExecutor, AsyncPipelineExecutor
from ..src.models import PipelineStep, UserRequest, TaskType


def test_pipeline_executor_initialization():
    """测试PipelineExecutor初始化"""
    config = {
        "timeout": 3600,
        "command_whitelist": []
    }
    executor = PipelineExecutor(config)
    assert executor is not None
    assert executor.timeout == 3600


def test_pipeline_executor_is_command_safe():
    """测试PipelineExecutor is_command_safe方法"""
    # 测试无白名单的情况
    config = {
        "timeout": 3600,
        "command_whitelist": []
    }
    executor = PipelineExecutor(config)
    assert executor._is_command_safe("echo hello") is True
    
    # 测试有白名单的情况
    config_with_whitelist = {
        "timeout": 3600,
        "command_whitelist": ["echo", "ls"]
    }
    executor_with_whitelist = PipelineExecutor(config_with_whitelist)
    assert executor_with_whitelist._is_command_safe("echo hello") is True
    assert executor_with_whitelist._is_command_safe("rm -rf /") is False


def test_pipeline_executor_fill_command(sample_pipeline_step):
    """测试PipelineExecutor fill_command方法"""
    config = {
        "timeout": 3600,
        "command_whitelist": []
    }
    executor = PipelineExecutor(config)
    
    context = {
        "input1": "sample1_1.fastq",
        "input2": "sample1_2.fastq",
        "output1": "sample1_1_trimmed.fastq",
        "output2": "sample1_2_trimmed.fastq",
        "adapter": "adapters.fa"
    }
    
    command = executor._fill_command(sample_pipeline_step, context)
    assert "sample1_1.fastq" in command
    assert "sample1_2.fastq" in command
    assert "sample1_1_trimmed.fastq" in command
    assert "sample1_2_trimmed.fastq" in command
    assert "adapters.fa" in command

@mock.patch('subprocess.run')
def test_pipeline_executor_run_step(mock_run, temp_dir, sample_pipeline_step):
    """测试PipelineExecutor run_step方法"""
    # 配置mock
    mock_run.return_value = mock.Mock(
        returncode=0,
        stdout=b"Command executed successfully",
        stderr=b""
    )
    
    config = {
        "timeout": 3600,
        "command_whitelist": []
    }
    executor = PipelineExecutor(config)
    
    context = {
        "input1": "sample1_1.fastq",
        "input2": "sample1_2.fastq",
        "output1": "sample1_1_trimmed.fastq",
        "output2": "sample1_2_trimmed.fastq",
        "adapter": "adapters.fa"
    }
    
    record = executor.run_step(sample_pipeline_step, context, temp_dir)
    assert record is not None
    assert record.status == "success"
    assert record.stdout == "Command executed successfully"
    assert mock_run.called

@mock.patch('subprocess.run')
def test_pipeline_executor_run_step_failure(mock_run, temp_dir, sample_pipeline_step):
    """测试PipelineExecutor run_step方法失败情况"""
    # 配置mock
    mock_run.return_value = mock.Mock(
        returncode=1,
        stdout=b"",
        stderr=b"Command failed"
    )
    
    config = {
        "timeout": 3600,
        "command_whitelist": []
    }
    executor = PipelineExecutor(config)
    
    context = {
        "input1": "sample1_1.fastq",
        "input2": "sample1_2.fastq",
        "output1": "sample1_1_trimmed.fastq",
        "output2": "sample1_2_trimmed.fastq",
        "adapter": "adapters.fa"
    }
    
    record = executor.run_step(sample_pipeline_step, context, temp_dir)
    assert record is not None
    assert record.status == "failed"
    assert record.stderr == "Command failed"

@pytest.mark.asyncio
@mock.patch('asyncio.create_subprocess_shell')
async def test_async_pipeline_executor_run_step(mock_create_subprocess, temp_dir, sample_pipeline_step):
    """测试AsyncPipelineExecutor run_step方法"""
    # 配置mock
    mock_process = mock.Mock()
    mock_process.communicate.return_value = (b"Command executed successfully", b"")
    mock_process.returncode = 0
    mock_create_subprocess.return_value = mock_process
    
    config = {
        "timeout": 3600,
        "command_whitelist": []
    }
    executor = AsyncPipelineExecutor(config)
    
    context = {
        "input1": "sample1_1.fastq",
        "input2": "sample1_2.fastq",
        "output1": "sample1_1_trimmed.fastq",
        "output2": "sample1_2_trimmed.fastq",
        "adapter": "adapters.fa"
    }
    
    record = await executor.run_step(sample_pipeline_step, context, temp_dir)
    assert record is not None
    assert record.status == "success"
    assert record.stdout == "Command executed successfully"
    assert mock_create_subprocess.called

@pytest.mark.asyncio
async def test_async_pipeline_executor_initialization():
    """测试AsyncPipelineExecutor初始化"""
    config = {
        "timeout": 3600,
        "command_whitelist": [],
        "max_concurrency": 4
    }
    executor = AsyncPipelineExecutor(config)
    assert executor is not None
    assert executor.timeout == 3600
    assert executor.max_concurrency == 4
