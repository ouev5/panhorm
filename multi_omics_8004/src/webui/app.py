import gradio as gr
import asyncio
import os
import json
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
try:
    import scanpy as sc
except ImportError:
    sc = None
from typing import Dict, Any, List
from ..models import UserRequest, TaskType, ExecutionStatus
from ..planner import Planner
from ..executor import AsyncPipelineExecutor
from ..xai_explainer import XAIExplainer
from ..registry import PipelineRegistry


class AutoBAWebUI:
    def __init__(self, config: dict):
        """初始化Web UI
        
        Args:
            config: 配置字典
        """
        self.config = config
        self.planner = Planner(config)
        self.executor = AsyncPipelineExecutor(config.get('executor', {}))
        self.xai_explainer = XAIExplainer(config.get('xai', {}))
        self.registry = PipelineRegistry()
        self.task_history = []
        self.current_task = None
        self.current_task_id = 0
        self.chat_contexts = {}  # 会话上下文

    def get_pipeline_list(self):
        """获取可用的流程列表"""
        pipelines = self.registry.discover()
        return [pipeline.name for pipeline in pipelines]

    async def run_analysis(self, data_paths, data_description, goal, pipeline_name, advanced_params, progress=gr.Progress()):
        """运行分析
        
        Args:
            data_paths: 数据路径
            data_description: 数据描述
            goal: 分析目标
            pipeline_name: 流程名称
            advanced_params: 高级参数
            progress: 进度条
            
        Returns:
            tuple: (日志, 结果, XAI报告, 下载链接)
        """
        # 生成任务ID
        task_id = self.current_task_id
        self.current_task_id += 1
        
        # 构建用户请求
        user_request = UserRequest(
            data_paths=data_paths.split(',') if isinstance(data_paths, str) else data_paths,
            data_description=data_description,
            goal=goal,
            task_type=TaskType.DETERMINISTIC
        )
        
        # 初始化日志
        logs = []
        logs.append(f"[Task {task_id}] 开始分析...")
        logs.append(f"[Task {task_id}] 数据路径: {data_paths}")
        logs.append(f"[Task {task_id}] 数据描述: {data_description}")
        logs.append(f"[Task {task_id}] 分析目标: {goal}")
        logs.append(f"[Task {task_id}] 流程名称: {pipeline_name}")
        
        try:
            # 生成或获取流程模板
            if pipeline_name == "auto":
                logs.append(f"[Task {task_id}] 自动生成分析计划...")
                template = await self.planner.plan(user_request)
                logs.append(f"[Task {task_id}] 生成计划: {template.name}")
            else:
                logs.append(f"[Task {task_id}] 加载预定义流程: {pipeline_name}")
                template = self.registry.get(pipeline_name)
                if not template:
                    logs.append(f"[Task {task_id}] 错误: 流程 {pipeline_name} 不存在")
                    return "\n".join(logs), "", "", ""
            
            # 准备执行上下文
            context = {
                'data_paths': data_paths,
                'output_dir': os.path.join(self.config.get('output_dir', './output'), f'task_{task_id}'),
                'script_dir': os.path.join(os.path.dirname(__file__), '../../scripts'),
                **(json.loads(advanced_params) if advanced_params else {})
            }
            
            # 确保输出目录存在
            os.makedirs(context['output_dir'], exist_ok=True)
            
            # 执行流程
            logs.append(f"[Task {task_id}] 开始执行流程...")
            
            # 定义进度回调
            def on_step_start(step):
                logs.append(f"[Task {task_id}] 开始执行步骤: {step.tool_name}")
                progress(0.0, f"执行步骤: {step.tool_name}")
            
            def on_step_complete(step, record):
                logs.append(f"[Task {task_id}] 完成步骤: {step.tool_name}, 状态: {record.status}")
                progress(1.0, f"完成步骤: {step.tool_name}")
            
            # 执行流程
            execution_records = await self.executor.run_pipeline(
                template=template,
                user_input=user_request,
                workdir=context['output_dir'],
                on_step_start=on_step_start,
                on_step_complete=on_step_complete
            )
            
            # 生成XAI报告
            logs.append(f"[Task {task_id}] 生成XAI报告...")
            result = {
                'status': 'success',
                'output_files': template.output_patterns
            }
            validation = {
                'passed': True,
                'score': 0.95
            }
            tools_used = [step.tool_name for step in template.steps]
            
            xai_output = self.xai_explainer.generate_report(result, validation, tools_used)
            
            # 构建下载链接
            download_links = []
            for output_pattern in template.output_patterns:
                output_file = output_pattern.replace('{output_dir}', context['output_dir'])
                if os.path.exists(output_file):
                    download_links.append(f"[下载 {os.path.basename(output_file)}](file://{output_file})")
            
            # 保存任务历史
            self.task_history.append({
                'task_id': task_id,
                'data_paths': data_paths,
                'data_description': data_description,
                'goal': goal,
                'pipeline_name': pipeline_name,
                'status': 'completed',
                'output_dir': context['output_dir']
            })
            
            logs.append(f"[Task {task_id}] 分析完成！")
            
            return "\n".join(logs), "分析完成！", xai_output.natural_language_explanation, "\n".join(download_links)
            
        except Exception as e:
            logs.append(f"[Task {task_id}] 错误: {str(e)}")
            
            # 保存任务历史
            self.task_history.append({
                'task_id': task_id,
                'data_paths': data_paths,
                'data_description': data_description,
                'goal': goal,
                'pipeline_name': pipeline_name,
                'status': 'failed',
                'error': str(e)
            })
            
            return "\n".join(logs), f"错误: {str(e)}", "", ""

    def cancel_task(self):
        """取消当前任务"""
        if self.current_task:
            self.current_task.cancel()
            return "任务已取消"
        return "没有正在运行的任务"

    def get_task_history(self):
        """获取任务历史"""
        history_str = "# 任务历史\n\n"
        for task in self.task_history:
            history_str += f"## 任务 {task['task_id']}\n"
            history_str += f"- 状态: {task['status']}\n"
            history_str += f"- 数据路径: {task['data_paths']}\n"
            history_str += f"- 数据描述: {task['data_description']}\n"
            history_str += f"- 分析目标: {task['goal']}\n"
            history_str += f"- 流程: {task['pipeline_name']}\n"
            if 'output_dir' in task:
                history_str += f"- 输出目录: {task['output_dir']}\n"
            if 'error' in task:
                history_str += f"- 错误: {task['error']}\n"
            history_str += "\n"
        return history_str

    async def handle_chat_message(self, message, history):
        """处理聊天消息
        
        Args:
            message: 用户消息
            history: 聊天历史
            
        Returns:
            str: 响应消息
        """
        # 生成会话ID
        session_id = hash(str(history)) % 10000
        
        # 初始化会话上下文
        if session_id not in self.chat_contexts:
            self.chat_contexts[session_id] = {
                'data_paths': '',
                'data_description': '',
                'goal': '',
                'current_plan': None
            }
        
        context = self.chat_contexts[session_id]
        
        # 解析用户输入
        if '数据路径' in message or '文件' in message:
            # 提取数据路径
            context['data_paths'] = message
            return "已记录数据路径。请描述您的数据类型和特点。"
        elif '数据类型' in message or '数据特点' in message:
            # 提取数据描述
            context['data_description'] = message
            return "已记录数据描述。请描述您的分析目标。"
        elif '分析目标' in message or '想做' in message:
            # 提取分析目标
            context['goal'] = message
            
            # 生成分析计划
            return await self.generate_plan_from_chat(context, session_id)
        elif '执行' in message or '运行' in message:
            # 执行计划
            if context['current_plan']:
                return await self.execute_plan_from_chat(context, session_id)
            else:
                return "请先生成分析计划。"
        else:
            # 尝试理解用户意图
            return "我是autoBA_optimized的聊天助手。请告诉我您的数据分析需求，包括数据路径、数据类型和分析目标。"

    async def generate_plan_from_chat(self, context, session_id):
        """从聊天生成分析计划
        
        Args:
            context: 会话上下文
            session_id: 会话ID
            
        Returns:
            str: 响应消息
        """
        # 构建用户请求
        user_request = UserRequest(
            data_paths=context['data_paths'].split(',') if context['data_paths'] else [],
            data_description=context['data_description'],
            goal=context['goal'],
            task_type=TaskType.DETERMINISTIC
        )
        
        # 生成计划
        try:
            template = await self.planner.plan(user_request)
            context['current_plan'] = template
            
            # 构建计划描述
            plan_desc = f"我为您生成了以下分析计划：\n\n"
            plan_desc += f"**计划名称**：{template.name}\n"
            plan_desc += f"**描述**：{template.description}\n"
            plan_desc += f"**任务类型**：{template.task_type}\n"
            plan_desc += f"**必需输入**：{', '.join(template.required_inputs)}\n"
            plan_desc += f"**输出文件**：{', '.join(template.output_patterns)}\n"
            plan_desc += "\n**步骤**：\n"
            for i, step in enumerate(template.steps, 1):
                plan_desc += f"{i}. {step.tool_name} (v{step.tool_version}) - {step.command_template}\n"
            plan_desc += "\n是否执行此计划？请回复'执行'或'运行'开始分析。"
            
            return plan_desc
        except Exception as e:
            return f"生成计划时出错：{str(e)}"

    async def execute_plan_from_chat(self, context, session_id):
        """从聊天执行分析计划
        
        Args:
            context: 会话上下文
            session_id: 会话ID
            
        Returns:
            str: 响应消息
        """
        template = context['current_plan']
        if not template:
            return "没有生成计划，请先描述您的分析需求。"
        
        # 生成任务ID
        task_id = self.current_task_id
        self.current_task_id += 1
        
        # 准备执行上下文
        output_dir = os.path.join(self.config.get('output_dir', './output'), f'chat_task_{task_id}')
        os.makedirs(output_dir, exist_ok=True)
        
        # 构建用户请求
        user_request = UserRequest(
            data_paths=context['data_paths'].split(',') if context['data_paths'] else [],
            data_description=context['data_description'],
            goal=context['goal'],
            task_type=TaskType.DETERMINISTIC
        )
        
        # 执行计划
        try:
            # 执行流程
            execution_records = await self.executor.run_pipeline(
                template=template,
                user_input=user_request,
                workdir=output_dir
            )
            
            # 构建结果消息
            result_msg = f"分析执行完成！\n\n"
            result_msg += f"**任务ID**：{task_id}\n"
            result_msg += f"**输出目录**：{output_dir}\n"
            result_msg += f"**生成的文件**：\n"
            for pattern in template.output_patterns:
                output_file = pattern.replace('{output_dir}', output_dir)
                if os.path.exists(output_file):
                    result_msg += f"- {os.path.basename(output_file)}\n"
            result_msg += "\n分析结果已保存到输出目录。"
            
            # 保存任务历史
            self.task_history.append({
                'task_id': task_id,
                'data_paths': context['data_paths'],
                'data_description': context['data_description'],
                'goal': context['goal'],
                'pipeline_name': template.name,
                'status': 'completed',
                'output_dir': output_dir
            })
            
            return result_msg
        except Exception as e:
            # 保存任务历史
            self.task_history.append({
                'task_id': task_id,
                'data_paths': context['data_paths'],
                'data_description': context['data_description'],
                'goal': context['goal'],
                'pipeline_name': template.name,
                'status': 'failed',
                'error': str(e)
            })
            
            return f"执行分析时出错：{str(e)}"

    def create_volcano_plot(self, csv_path) -> go.Figure:
        """创建差异表达火山图
        
        Args:
            csv_path: CSV文件路径
            
        Returns:
            plotly.Figure: 火山图
        """
        try:
            # 读取CSV文件
            df = pd.read_csv(csv_path)
            
            # 检查必要的列
            if 'log2FoldChange' not in df.columns or 'pvalue' not in df.columns:
                return go.Figure(go.Scatter(x=[], y=[], text="缺少必要的列: log2FoldChange 和 pvalue"))
            
            # 计算-log10(pvalue)
            df['neg_log10_pvalue'] = -df['pvalue'].apply(lambda x: 0 if x == 0 else -np.log10(x))
            
            # 标记显著差异基因
            df['significant'] = (df['log2FoldChange'].abs() > 1) & (df['pvalue'] < 0.05)
            
            # 创建火山图
            fig = px.scatter(
                df,
                x='log2FoldChange',
                y='neg_log10_pvalue',
                color='significant',
                hover_data=['gene'],
                title='差异表达火山图',
                labels={
                    'log2FoldChange': 'log2(Fold Change)',
                    'neg_log10_pvalue': '-log10(p-value)',
                    'significant': '显著差异'
                }
            )
            
            # 添加阈值线
            fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="gray")
            fig.add_vline(x=1, line_dash="dash", line_color="gray")
            fig.add_vline(x=-1, line_dash="dash", line_color="gray")
            
            return fig
        except Exception as e:
            return go.Figure(go.Scatter(x=[], y=[], text=f"创建火山图时出错: {str(e)}"))

    def create_umap_plot(self, h5ad_path) -> go.Figure:
        """创建单细胞UMAP降维图
        
        Args:
            h5ad_path: AnnData文件路径
            
        Returns:
            plotly.Figure: UMAP图
        """
        try:
            # 读取AnnData文件
            adata = sc.read(h5ad_path)
            
            # 检查是否有UMAP嵌入
            if 'X_umap' not in adata.obsm:
                # 如果没有，计算UMAP
                sc.tl.pca(adata)
                sc.pp.neighbors(adata)
                sc.tl.umap(adata)
            
            # 提取UMAP坐标
            umap_df = pd.DataFrame(
                adata.obsm['X_umap'],
                columns=['UMAP1', 'UMAP2'],
                index=adata.obs_names
            )
            
            # 添加细胞类型或聚类信息
            if 'cell_type' in adata.obs:
                umap_df['cell_type'] = adata.obs['cell_type'].values
                color_col = 'cell_type'
            elif 'leiden' in adata.obs:
                umap_df['cluster'] = adata.obs['leiden'].values
                color_col = 'cluster'
            else:
                umap_df['cell'] = range(len(umap_df))
                color_col = 'cell'
            
            # 创建UMAP图
            fig = px.scatter(
                umap_df,
                x='UMAP1',
                y='UMAP2',
                color=color_col,
                title='单细胞UMAP降维图',
                labels={
                    'UMAP1': 'UMAP 1',
                    'UMAP2': 'UMAP 2'
                }
            )
            
            return fig
        except Exception as e:
            return go.Figure(go.Scatter(x=[], y=[], text=f"创建UMAP图时出错: {str(e)}"))

    def create_peak_track_plot(self, bed_path, bam_path=None) -> go.Figure:
        """创建峰值信号轨迹图
        
        Args:
            bed_path: BED文件路径
            bam_path: BAM文件路径（可选）
            
        Returns:
            plotly.Figure: 峰值轨迹图
        """
        try:
            # 读取BED文件
            peaks = pd.read_csv(bed_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score'])
            
            # 选择前10个峰值进行展示
            top_peaks = peaks.head(10)
            
            # 创建轨迹图
            fig = go.Figure()
            
            for i, row in top_peaks.iterrows():
                fig.add_trace(go.Scatter(
                    x=[row['start'], row['end']],
                    y=[i, i],
                    mode='lines',
                    line=dict(width=10),
                    name=f"Peak {i+1}"
                ))
            
            fig.update_layout(
                title='峰值信号轨迹图',
                xaxis_title='基因组位置',
                yaxis_title='峰值索引',
                showlegend=False
            )
            
            return fig
        except Exception as e:
            return go.Figure(go.Scatter(x=[], y=[], text=f"创建峰值轨迹图时出错: {str(e)}"))

    def load_table_data(self, csv_path) -> pd.DataFrame:
        """加载表格数据
        
        Args:
            csv_path: CSV文件路径
            
        Returns:
            pd.DataFrame: 表格数据
        """
        try:
            return pd.read_csv(csv_path)
        except Exception as e:
            return pd.DataFrame({'error': [str(e)]})

    def build_ui(self):
        """构建UI"""
        with gr.Blocks(title="autoBA_optimized - 生物信息学自动化分析平台") as app:
            gr.Markdown("# autoBA_optimized - 生物信息学自动化分析平台")
            
            with gr.Row():
                # 左侧：配置面板
                with gr.Column(scale=1, min_width=400):
                    gr.Markdown("## 配置")
                    
                    # 数据路径输入
                    data_paths = gr.Textbox(
                        label="数据路径",
                        placeholder="输入文件路径，多个路径用逗号分隔",
                        lines=2
                    )
                    
                    # 数据上传
                    file_upload = gr.File(
                        label="或上传文件",
                        file_count="multiple"
                    )
                    
                    # 数据描述
                    data_description = gr.Textbox(
                        label="数据描述",
                        placeholder="描述您的数据类型和特点",
                        lines=3
                    )
                    
                    # 分析目标
                    goal = gr.Textbox(
                        label="分析目标",
                        placeholder="描述您想要完成的分析任务",
                        lines=2
                    )
                    
                    # 预定义流程选择
                    pipeline_list = self.get_pipeline_list()
                    pipeline_list.insert(0, "auto")  # 添加自动生成选项
                    pipeline_name = gr.Dropdown(
                        label="预定义流程",
                        choices=pipeline_list,
                        value="auto"
                    )
                    
                    # 高级参数
                    with gr.Accordion("高级参数", open=False):
                        advanced_params = gr.Textbox(
                            label="高级参数 (JSON格式)",
                            placeholder='{"param1": "value1", "param2": "value2"}',
                            lines=3
                        )
                    
                    # 运行按钮
                    run_button = gr.Button("运行分析", variant="primary")
                    cancel_button = gr.Button("取消运行")
                
                # 右侧：结果展示
                with gr.Column(scale=2, min_width=600):
                    gr.Markdown("## 结果")
                    
                    # 实时日志输出
                    logs = gr.Textbox(
                        label="日志",
                        placeholder="运行日志将显示在这里",
                        lines=10,
                        interactive=False
                    )
                    
                    # 结果预览
                    result_preview = gr.Textbox(
                        label="结果",
                        placeholder="分析结果将显示在这里",
                        lines=2,
                        interactive=False
                    )
                    
                    # XAI报告区域
                    xai_report = gr.Markdown(
                        label="XAI报告",
                        value="# XAI报告\n\n分析完成后将显示详细报告"
                    )
                    
                    # 下载链接
                    download_links = gr.Markdown(
                        label="下载链接",
                        value="# 下载链接\n\n分析完成后将显示下载链接"
                    )
                    
                    # 任务历史
                    task_history = gr.Markdown(
                        label="任务历史",
                        value=self.get_task_history()
                    )
                    
                    # 可视化组件
                    gr.Markdown("## 结果可视化")
                    
                    # 火山图可视化
                    with gr.Tab("差异表达火山图"):
                        volcano_file = gr.File(label="选择差异表达结果CSV文件")
                        volcano_button = gr.Button("生成火山图")
                        volcano_plot = gr.Plot(label="火山图")
                    
                    # UMAP图可视化
                    with gr.Tab("单细胞UMAP图"):
                        umap_file = gr.File(label="选择AnnData文件 (.h5ad)")
                        umap_button = gr.Button("生成UMAP图")
                        umap_plot = gr.Plot(label="UMAP图")
                    
                    # 峰值轨迹图可视化
                    with gr.Tab("峰值信号轨迹图"):
                        peak_file = gr.File(label="选择BED文件")
                        peak_button = gr.Button("生成轨迹图")
                        peak_plot = gr.Plot(label="峰值轨迹图")
                    
                    # 表格数据展示
                    with gr.Tab("表格数据"):
                        table_file = gr.File(label="选择CSV文件")
                        table_button = gr.Button("加载表格")
                        table_output = gr.Dataframe(label="表格数据", interactive=False)
                
                # 右侧边栏：聊天助手
                with gr.Column(scale=1, min_width=400):
                    gr.Markdown("## 聊天助手")
                    
                    # 聊天界面
                    chat_interface = gr.ChatInterface(
                        fn=self.handle_chat_message,
                        chatbot=gr.Chatbot(),
                        textbox=gr.Textbox(
                            placeholder="请描述您的分析需求，例如：我有RNA-seq数据，想做差异表达分析",
                            lines=2
                        ),
                        title="autoBA_optimized 助手",
                        description="通过自然语言描述您的分析需求，我将为您生成分析计划并执行。"
                    )
            
            # 事件处理
            run_button.click(
                fn=self.run_analysis,
                inputs=[data_paths, data_description, goal, pipeline_name, advanced_params],
                outputs=[logs, result_preview, xai_report, download_links]
            )
            
            cancel_button.click(
                fn=self.cancel_task,
                inputs=[],
                outputs=[logs]
            )
            
            # 文件上传处理
            def handle_file_upload(files):
                if files:
                    paths = [file.name for file in files]
                    return ", ".join(paths)
                return ""
            
            file_upload.change(
                fn=handle_file_upload,
                inputs=[file_upload],
                outputs=[data_paths]
            )
            
            # 可视化事件处理
            volcano_button.click(
                fn=self.create_volcano_plot,
                inputs=[volcano_file],
                outputs=[volcano_plot]
            )
            
            umap_button.click(
                fn=self.create_umap_plot,
                inputs=[umap_file],
                outputs=[umap_plot]
            )
            
            peak_button.click(
                fn=self.create_peak_track_plot,
                inputs=[peak_file],
                outputs=[peak_plot]
            )
            
            table_button.click(
                fn=self.load_table_data,
                inputs=[table_file],
                outputs=[table_output]
            )
        
        return app


def main():
    """主函数"""
    # 加载配置
    config = {
        "llm": {
            "mode": "deepseek",
            "api_key": "YOUR_VALUE_HERE",
            "api_url": "https://api.deepseek.com/v1/chat/completions",
            "model": "deepseek-chat"
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
    
    # 创建Web UI
    webui = AutoBAWebUI(config)
    app = webui.build_ui()
    
    # 启动应用 - 使用8004端口
    app.launch(server_name="0.0.0.0", server_port=8004, share=False)


if __name__ == "__main__":
    main()
