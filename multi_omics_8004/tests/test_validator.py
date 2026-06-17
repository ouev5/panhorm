import pytest
import os
import tempfile
import pandas as pd
from ..src.validator import RuleValidator, LLMValidator


def test_rule_validator_check_file_not_empty(temp_dir):
    """测试RuleValidator check_file_not_empty方法"""
    validator = RuleValidator({})
    
    # 测试空文件
    empty_file = os.path.join(temp_dir, "empty.txt")
    with open(empty_file, "w") as f:
        pass
    result = validator.check_file_not_empty(empty_file)
    assert not result.passed
    assert result.score < 0.5
    
    # 测试非空文件
    non_empty_file = os.path.join(temp_dir, "non_empty.txt")
    with open(non_empty_file, "w") as f:
        f.write("Hello, world!")
    result = validator.check_file_not_empty(non_empty_file)
    assert result.passed
    assert result.score >= 0.9
    
    # 测试不存在的文件
    non_existent_file = os.path.join(temp_dir, "non_existent.txt")
    result = validator.check_file_not_empty(non_existent_file)
    assert not result.passed
    assert result.score < 0.5


def test_rule_validator_check_csv_columns(temp_dir):
    """测试RuleValidator check_csv_columns方法"""
    validator = RuleValidator({})
    
    # 测试包含所有必需列的CSV
    csv_file = os.path.join(temp_dir, "test.csv")
    df = pd.DataFrame({
        "gene": ["A", "B", "C"],
        "log2FoldChange": [1.0, 2.0, -1.0],
        "pvalue": [0.01, 0.001, 0.05]
    })
    df.to_csv(csv_file, index=False)
    
    required_columns = ["gene", "log2FoldChange", "pvalue"]
    result = validator.check_csv_columns(csv_file, required_columns)
    assert result.passed
    assert result.score >= 0.9
    
    # 测试缺少列的CSV
    df_missing = pd.DataFrame({
        "gene": ["A", "B", "C"],
        "log2FoldChange": [1.0, 2.0, -1.0]
    })
    df_missing.to_csv(csv_file, index=False)
    result = validator.check_csv_columns(csv_file, required_columns)
    assert not result.passed
    assert result.score < 0.5


def test_rule_validator_check_log2fc_distribution(temp_dir):
    """测试RuleValidator check_log2fc_distribution方法"""
    validator = RuleValidator({})
    
    # 测试正常分布的log2FC
    csv_file = os.path.join(temp_dir, "test.csv")
    df = pd.DataFrame({
        "gene": [f"gene{i}" for i in range(100)],
        "log2FoldChange": [i/10 for i in range(-50, 50)],
        "pvalue": [0.01 for _ in range(100)]
    })
    df.to_csv(csv_file, index=False)
    
    result = validator.check_log2fc_distribution(csv_file)
    assert result.passed
    assert result.score >= 0.8


def test_rule_validator_check_bam_index(temp_dir):
    """测试RuleValidator check_bam_index方法"""
    validator = RuleValidator({})
    
    # 测试有索引的BAM文件
    bam_file = os.path.join(temp_dir, "test.bam")
    bai_file = os.path.join(temp_dir, "test.bam.bai")
    
    # 创建空文件
    with open(bam_file, "w") as f:
        pass
    with open(bai_file, "w") as f:
        pass
    
    result = validator.check_bam_index(bam_file)
    assert result.passed
    assert result.score >= 0.9
    
    # 测试没有索引的BAM文件
    os.remove(bai_file)
    result = validator.check_bam_index(bam_file)
    assert not result.passed
    assert result.score < 0.5


def test_rule_validator_check_peak_count(temp_dir):
    """测试RuleValidator check_peak_count方法"""
    validator = RuleValidator({})
    
    # 测试峰值数量足够的BED文件
    bed_file = os.path.join(temp_dir, "peaks.bed")
    with open(bed_file, "w") as f:
        for i in range(150):
            f.write(f"chr1\t{i*100}\t{i*100+100}\tpeak{i}\t100\n")
    
    result = validator.check_peak_count(bed_file)
    assert result.passed
    assert result.score >= 0.9
    
    # 测试峰值数量不足的BED文件
    bed_file_small = os.path.join(temp_dir, "peaks_small.bed")
    with open(bed_file_small, "w") as f:
        for i in range(50):
            f.write(f"chr1\t{i*100}\t{i*100+100}\tpeak{i}\t100\n")
    
    result = validator.check_peak_count(bed_file_small)
    assert not result.passed
    assert result.score < 0.5


def test_llm_validator_initialization():
    """测试LLMValidator初始化"""
    config = {
        "llm": {
            "mode": "mock"
        }
    }
    validator = LLMValidator(config)
    assert validator is not None


def test_llm_validator_generate_prompt():
    """测试LLMValidator生成prompt"""
    config = {
        "llm": {
            "mode": "mock"
        }
    }
    validator = LLMValidator(config)
    
    # 测试生成prompt方法（内部方法）
    task_name = "RNA-seq差异表达"
    content = "基因A的log2FC为2.5，p值为0.01"
    
    # 这里我们测试prompt的构造逻辑
    # 由于prompt生成是内部方法，我们通过间接方式测试
    assert task_name in validator._construct_prompt(task_name, content)
    assert content in validator._construct_prompt(task_name, content)


def test_llm_validator_parse_response():
    """测试LLMValidator解析响应"""
    config = {
        "llm": {
            "mode": "mock"
        }
    }
    validator = LLMValidator(config)
    
    # 测试解析有效响应
    valid_response = '{"passed": true, "score": 0.9, "issues": [], "suggestion": "分析结果正常"}'
    result = validator._parse_llm_response(valid_response)
    assert result is not None
    assert result["passed"] is True
    assert result["score"] == 0.9
    
    # 测试解析无效响应
    invalid_response = "这不是一个有效的JSON"
    result = validator._parse_llm_response(invalid_response)
    assert result is None
