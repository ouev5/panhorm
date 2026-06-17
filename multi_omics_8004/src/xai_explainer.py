import json
from typing import Dict, Any, List, Optional
from models import XAIOutput
from llm_client import LLMPlannerClient
from rag_retriever import RAGRetriever


class XAIExplainer:
    def __init__(self, config: dict):
        """初始化XAI解释器
        
        Args:
            config: 配置字典
        """
        self.config = config
        self.llm_client = LLMPlannerClient(config.get('llm', {}))
        self.rag_retriever = RAGRetriever(config.get('rag_persist_dir', './rag_index'))
        self.language = config.get('language', 'zh')  # zh or en

    def explain_decision(self, decision: Dict[str, Any], context: Dict[str, Any]) -> str:
        """解释为什么选择某个流程
        
        Args:
            decision: 决策数据
            context: 上下文数据
            
        Returns:
            str: 解释文本
        """
        # 构建Prompt
        if self.language == 'zh':
            system_message = {
                "role": "system",
                "content": "你是一个生物信息学分析专家，负责解释分析流程的选择决策。请根据提供的决策数据和上下文，生成详细的解释，说明为什么选择了该流程。"
            }
            user_message = {
                "role": "user",
                "content": f"决策数据：\n{json.dumps(decision, ensure_ascii=False, indent=2)}\n\n上下文数据：\n{json.dumps(context, ensure_ascii=False, indent=2)}\n\n请生成一个详细的解释，说明为什么选择了该流程，包括考虑的因素和证据。"
            }
        else:
            system_message = {
                "role": "system",
                "content": "You are a bioinformatics analysis expert responsible for explaining the selection decision of an analysis pipeline. Please generate a detailed explanation based on the provided decision data and context, explaining why this pipeline was chosen."
            }
            user_message = {
                "role": "user",
                "content": f"Decision data: \n{json.dumps(decision, indent=2)}\n\nContext data: \n{json.dumps(context, indent=2)}\n\nPlease generate a detailed explanation of why this pipeline was chosen, including the factors considered and evidence."
            }
        
        # 调用LLM
        response = self.llm_client.generate([system_message, user_message])
        
        # 提取内容
        if response and 'choices' in response:
            return response['choices'][0]['message']['content']
        return "无法生成解释"

    def explain_error(self, error_log: str, suggestion: str) -> str:
        """解释错误原因
        
        Args:
            error_log: 错误日志
            suggestion: 建议
            
        Returns:
            str: 解释文本
        """
        # 构建Prompt
        if self.language == 'zh':
            system_message = {
                "role": "system",
                "content": "你是一个生物信息学分析专家，负责解释分析过程中的错误原因。请根据提供的错误日志和建议，生成详细的解释，说明错误的原因和如何解决。"
            }
            user_message = {
                "role": "user",
                "content": f"错误日志：\n{error_log}\n\n建议：\n{suggestion}\n\n请生成一个详细的解释，说明错误的原因和如何解决。"
            }
        else:
            system_message = {
                "role": "system",
                "content": "You are a bioinformatics analysis expert responsible for explaining error causes in analysis processes. Please generate a detailed explanation based on the provided error log and suggestion, explaining the error cause and how to resolve it."
            }
            user_message = {
                "role": "user",
                "content": f"Error log: \n{error_log}\n\nSuggestion: \n{suggestion}\n\nPlease generate a detailed explanation of the error cause and how to resolve it."
            }
        
        # 调用LLM
        response = self.llm_client.generate([system_message, user_message])
        
        # 提取内容
        if response and 'choices' in response:
            return response['choices'][0]['message']['content']
        return "无法生成错误解释"

    def generate_report(self, result: Dict[str, Any], validation: Dict[str, Any], tools_used: List[str]) -> XAIOutput:
        """生成XAI报告
        
        Args:
            result: 分析结果
            validation: 验证结果
            tools_used: 使用的工具列表
            
        Returns:
            XAIOutput: XAI输出
        """
        # 检索相关文档作为证据
        query = f"分析结果：{result.get('status', '')}，使用工具：{', '.join(tools_used)}"
        evidence_sources = self.rag_retriever.retrieve(query, top_k=3)
        
        # 构建证据来源列表
        evidence_list = []
        for i, doc in enumerate(evidence_sources):
            evidence_list.append({
                'id': i + 1,
                'source': doc.metadata.get('source', 'Unknown'),
                'content': doc.page_content[:300] + '...' if len(doc.page_content) > 300 else doc.page_content
            })
        
        # 构建Prompt
        if self.language == 'zh':
            system_message = {
                "role": "system",
                "content": "你是一个生物信息学分析专家，负责生成分析结果的解释报告。请根据提供的分析结果、验证结果和使用的工具，生成一个详细的Markdown格式报告，包括分析过程、结果解释、证据来源和建议。"
            }
            user_message = {
                "role": "user",
                "content": f"分析结果：\n{json.dumps(result, ensure_ascii=False, indent=2)}\n\n验证结果：\n{json.dumps(validation, ensure_ascii=False, indent=2)}\n\n使用的工具：\n{', '.join(tools_used)}\n\n证据来源：\n{json.dumps(evidence_list, ensure_ascii=False, indent=2)}\n\n请生成一个详细的Markdown格式报告，包括分析过程、结果解释、证据来源和建议。"
            }
        else:
            system_message = {
                "role": "system",
                "content": "You are a bioinformatics analysis expert responsible for generating explanation reports for analysis results. Please generate a detailed Markdown format report based on the provided analysis results, validation results, and tools used, including analysis process, result explanation, evidence sources, and suggestions."
            }
            user_message = {
                "role": "user",
                "content": f"Analysis results: \n{json.dumps(result, indent=2)}\n\nValidation results: \n{json.dumps(validation, indent=2)}\n\nTools used: \n{', '.join(tools_used)}\n\nEvidence sources: \n{json.dumps(evidence_list, indent=2)}\n\nPlease generate a detailed Markdown format report, including analysis process, result explanation, evidence sources, and suggestions."
            }
        
        # 调用LLM
        response = self.llm_client.generate([system_message, user_message])
        
        # 提取内容
        natural_language_explanation = "无法生成解释报告"
        if response and 'choices' in response:
            natural_language_explanation = response['choices'][0]['message']['content']
        
        # 计算置信度
        confidence = 0.85  # 示例值，实际应根据验证结果计算
        
        # 生成特征贡献（示例）
        feature_contributions = []
        for tool in tools_used:
            feature_contributions.append({
                'feature': tool,
                'contribution': 1.0 / len(tools_used)
            })
        
        # 创建XAIOutput
        xai_output = XAIOutput(
            confidence=confidence,
            evidence_sources=evidence_list,
            feature_contributions=feature_contributions,
            natural_language_explanation=natural_language_explanation
        )
        
        return xai_output

    def visualize_feature_importance(self, shap_values: Dict[str, Any]) -> str:
        """SHAP可视化
        
        Args:
            shap_values: SHAP值
            
        Returns:
            str: 可视化结果（HTML或Markdown）
        """
        # 这里简化实现，实际应使用SHAP库生成可视化
        if self.language == 'zh':
            visualization = "## 特征重要性可视化\n\n"
            visualization += "### SHAP值分析\n\n"
            for feature, value in shap_values.items():
                visualization += f"- **{feature}**: {value}\n"
            visualization += "\n### 解释\n\n特征重要性反映了每个特征对模型预测的贡献程度。"
        else:
            visualization = "## Feature Importance Visualization\n\n"
            visualization += "### SHAP Value Analysis\n\n"
            for feature, value in shap_values.items():
                visualization += f"- **{feature}**: {value}\n"
            visualization += "\n### Explanation\n\nFeature importance reflects the contribution of each feature to the model's prediction."
        
        return visualization
