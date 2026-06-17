import subprocess
import os
import sys
import time
import asyncio
from datetime import datetime
from typing import Dict, Any, List, Optional, Callable
from models import PipelineStep, PipelineTemplate, ExecutionRecord, ExecutionStatus


class PipelineExecutor:
    def __init__(self, config: dict):
        """初始化执行器
        
        Args:
            config: 配置字典，包含超时、白名单命令等
        """
        self.config = config
        self.timeout = config.get('timeout', 3600)  # 默认超时1小时
        self.command_whitelist = config.get('command_whitelist', [])
        self.is_windows = sys.platform.startswith('win32')

    def _is_command_safe(self, command: str) -> bool:
        """检查命令是否安全
        
        Args:
            command: 要检查的命令
            
        Returns:
            bool: 命令是否安全
        """
        # 如果没有配置白名单，则默认允许所有命令
        if not self.command_whitelist:
            return True
        
        # 检查命令是否在白名单中
        for allowed_command in self.command_whitelist:
            if command.startswith(allowed_command):
                return True
        
        return False

    def _fill_command(self, step: PipelineStep, context: dict) -> str:
        """渲染命令模板
        
        Args:
            step: 管道步骤
            context: 上下文变量
            
        Returns:
            str: 渲染后的命令
        """
        return step.render_command(context)

    def run_step(self, step: PipelineStep, context: dict, workdir: str) -> ExecutionRecord:
        """同步执行单步
        
        Args:
            step: 要执行的步骤
            context: 上下文变量
            workdir: 工作目录
            
        Returns:
            ExecutionRecord: 执行记录
        """
        execution_id = f"step_{step.step_id}_{int(time.time())}"
        record = ExecutionRecord(
            execution_id=execution_id,
            pipeline_name=f"step_{step.step_id}_{step.tool_name}",
            status=ExecutionStatus.RUNNING
        )
        
        try:
            # 渲染命令
            command = self._fill_command(step, context)
            
            # 安全检查
            if not self._is_command_safe(command):
                record.status = ExecutionStatus.FAILED
                record.error_message = f"命令不在白名单中: {command}"
                record.end_time = datetime.now()
                return record
            
            # 执行命令
            start_time = datetime.now()
            
            # 处理Windows和Linux的兼容性
            if self.is_windows:
                # Windows使用shell=True
                result = subprocess.run(
                    command,
                    shell=True,
                    cwd=workdir,
                    capture_output=True,
                    text=True,
                    timeout=self.timeout
                )
            else:
                # Linux使用shell=False，需要拆分命令
                import shlex
                cmd_args = shlex.split(command)
                result = subprocess.run(
                    cmd_args,
                    shell=False,
                    cwd=workdir,
                    capture_output=True,
                    text=True,
                    timeout=self.timeout
                )
            
            end_time = datetime.now()
            
            # 记录执行结果
            step_execution = {
                'step_id': step.step_id,
                'tool_name': step.tool_name,
                'command': command,
                'stdout': result.stdout,
                'stderr': result.stderr,
                'returncode': result.returncode,
                'start_time': start_time.isoformat(),
                'end_time': end_time.isoformat()
            }
            
            record.step_executions.append(step_execution)
            record.end_time = end_time
            
            if result.returncode == 0:
                record.status = ExecutionStatus.SUCCESS
            else:
                record.status = ExecutionStatus.FAILED
                record.error_message = f"命令执行失败: {result.stderr}"
                
        except subprocess.TimeoutExpired:
            record.status = ExecutionStatus.FAILED
            record.error_message = f"命令执行超时（{self.timeout}秒）"
            record.end_time = datetime.now()
        except Exception as e:
            record.status = ExecutionStatus.FAILED
            record.error_message = f"执行错误: {str(e)}"
            record.end_time = datetime.now()
        
        return record

    def run_pipeline(self, template: PipelineTemplate, user_input: dict, workdir: str) -> List[ExecutionRecord]:
        """执行整个流程
        
        Args:
            template: 管道模板
            user_input: 用户输入
            workdir: 工作目录
            
        Returns:
            List[ExecutionRecord]: 执行记录列表
        """
        execution_records = []
        context = user_input.copy()
        
        # 确保工作目录存在
        os.makedirs(workdir, exist_ok=True)
        
        # 执行每个步骤
        for step in template.steps:
            # 更新上下文，添加前一步的输出作为输入
            if execution_records:
                last_record = execution_records[-1]
                if last_record.status == ExecutionStatus.SUCCESS:
                    # 可以在这里添加逻辑，将前一步的输出作为下一步的输入
                    pass
                else:
                    # 如果前一步失败，停止执行
                    break
            
            # 执行当前步骤
            record = self.run_step(step, context, workdir)
            execution_records.append(record)
            
            # 如果步骤失败，停止执行
            if record.status != ExecutionStatus.SUCCESS:
                break
        
        return execution_records


class AsyncPipelineExecutor:
    def __init__(self, config: dict):
        """初始化异步执行器
        
        Args:
            config: 配置字典，包含超时、白名单命令等
        """
        self.config = config
        self.timeout = config.get('timeout', 3600)  # 默认超时1小时
        self.command_whitelist = config.get('command_whitelist', [])
        self.max_concurrency = config.get('max_concurrency', 4)  # 默认最大并发数4
        self.semaphore = asyncio.Semaphore(self.max_concurrency)
        self.is_windows = sys.platform.startswith('win32')

    def _is_command_safe(self, command: str) -> bool:
        """检查命令是否安全
        
        Args:
            command: 要检查的命令
            
        Returns:
            bool: 命令是否安全
        """
        # 如果没有配置白名单，则默认允许所有命令
        if not self.command_whitelist:
            return True
        
        # 检查命令是否在白名单中
        for allowed_command in self.command_whitelist:
            if command.startswith(allowed_command):
                return True
        
        return False

    def _fill_command(self, step: PipelineStep, context: dict) -> str:
        """渲染命令模板
        
        Args:
            step: 管道步骤
            context: 上下文变量
            
        Returns:
            str: 渲染后的命令
        """
        return step.render_command(context)

    async def run_step(self, step: PipelineStep, context: dict, workdir: str, 
                      on_step_start: Optional[Callable] = None, 
                      on_step_complete: Optional[Callable] = None) -> ExecutionRecord:
        """异步执行单步
        
        Args:
            step: 要执行的步骤
            context: 上下文变量
            workdir: 工作目录
            on_step_start: 步骤开始回调
            on_step_complete: 步骤完成回调
            
        Returns:
            ExecutionRecord: 执行记录
        """
        execution_id = f"step_{step.step_id}_{int(time.time())}"
        record = ExecutionRecord(
            execution_id=execution_id,
            pipeline_name=f"step_{step.step_id}_{step.tool_name}",
            status=ExecutionStatus.RUNNING
        )
        
        # 调用开始回调
        if on_step_start:
            await on_step_start(step, record)
        
        try:
            # 渲染命令
            command = self._fill_command(step, context)
            
            # 安全检查
            if not self._is_command_safe(command):
                record.status = ExecutionStatus.FAILED
                record.error_message = f"命令不在白名单中: {command}"
                record.end_time = datetime.now()
                if on_step_complete:
                    await on_step_complete(step, record)
                return record
            
            # 执行命令
            start_time = datetime.now()
            
            # 使用Semaphore控制并发
            async with self.semaphore:
                # 使用asyncio.timeout上下文管理器
                async with asyncio.timeout(self.timeout):
                    # 使用asyncio.create_subprocess_shell
                    process = await asyncio.create_subprocess_shell(
                        command,
                        cwd=workdir,
                        stdout=asyncio.subprocess.PIPE,
                        stderr=asyncio.subprocess.PIPE,
                        text=True
                    )
                    
                    # 等待进程完成
                    stdout, stderr = await process.communicate()
            
            end_time = datetime.now()
            
            # 记录执行结果
            step_execution = {
                'step_id': step.step_id,
                'tool_name': step.tool_name,
                'command': command,
                'stdout': stdout,
                'stderr': stderr,
                'returncode': process.returncode,
                'start_time': start_time.isoformat(),
                'end_time': end_time.isoformat()
            }
            
            record.step_executions.append(step_execution)
            record.end_time = end_time
            
            if process.returncode == 0:
                record.status = ExecutionStatus.SUCCESS
            else:
                record.status = ExecutionStatus.FAILED
                record.error_message = f"命令执行失败: {stderr}"
                
        except asyncio.TimeoutError:
            record.status = ExecutionStatus.FAILED
            record.error_message = f"命令执行超时（{self.timeout}秒）"
            record.end_time = datetime.now()
        except asyncio.CancelledError:
            record.status = ExecutionStatus.ABORTED
            record.error_message = "任务被取消"
            record.end_time = datetime.now()
        except Exception as e:
            record.status = ExecutionStatus.FAILED
            record.error_message = f"执行错误: {str(e)}"
            record.end_time = datetime.now()
        
        # 调用完成回调
        if on_step_complete:
            await on_step_complete(step, record)
        
        return record

    async def run_pipeline(self, template: PipelineTemplate, user_input: dict, workdir: str, 
                         on_step_start: Optional[Callable] = None, 
                         on_step_complete: Optional[Callable] = None) -> List[ExecutionRecord]:
        """异步执行整个流程
        
        Args:
            template: 管道模板
            user_input: 用户输入
            workdir: 工作目录
            on_step_start: 步骤开始回调
            on_step_complete: 步骤完成回调
            
        Returns:
            List[ExecutionRecord]: 执行记录列表
        """
        execution_records = []
        context = user_input.copy()
        
        # 确保工作目录存在
        os.makedirs(workdir, exist_ok=True)
        
        # 使用asyncio.taskgroup管理任务组
        async with asyncio.TaskGroup() as tg:
            tasks = []
            
            # 为每个步骤创建任务
            for step in template.steps:
                # 执行当前步骤
                task = tg.create_task(
                    self.run_step(step, context, workdir, on_step_start, on_step_complete)
                )
                tasks.append(task)
            
            # 等待所有任务完成
            for task in tasks:
                try:
                    record = await task
                    execution_records.append(record)
                except asyncio.CancelledError:
                    # 任务被取消
                    pass
        
        return execution_records

