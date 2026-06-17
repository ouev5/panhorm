import pytest
from ..src.arbiter import Arbiter, AdaptiveIterator
from ..src.models import ValidationResult


def test_arbiter_initialization():
    """测试Arbiter初始化"""
    config = {
        "max_retries": 3,
        "score_improvement_threshold": 0.1
    }
    arbiter = Arbiter(config)
    assert arbiter is not None
    assert arbiter.max_retries == 3
    assert arbiter.score_improvement_threshold == 0.1


def test_arbiter_decide_next_action_success():
    """测试Arbiter decide_next_action方法 - 成功情况"""
    config = {
        "max_retries": 3,
        "score_improvement_threshold": 0.1
    }
    arbiter = Arbiter(config)
    
    # 验证通过的情况
    history = []
    validation = ValidationResult(
        passed=True,
        score=0.9,
        issues=[],
        suggestion="验证通过"
    )
    
    decision = arbiter.decide_next_action(history, validation)
    assert decision == "success"


def test_arbiter_decide_next_action_abort():
    """测试Arbiter decide_next_action方法 - 中止情况"""
    config = {
        "max_retries": 3,
        "score_improvement_threshold": 0.1
    }
    arbiter = Arbiter(config)
    
    # 重试次数达到最大值的情况
    history = [
        {"status": "failed", "score": 0.5},
        {"status": "failed", "score": 0.6},
        {"status": "failed", "score": 0.7}
    ]
    validation = ValidationResult(
        passed=False,
        score=0.75,
        issues=["测试错误"],
        suggestion="修改参数"
    )
    
    decision = arbiter.decide_next_action(history, validation)
    assert decision == "abort"


def test_arbiter_decide_next_action_fallback():
    """测试Arbiter decide_next_action方法 - 回退情况"""
    config = {
        "max_retries": 3,
        "score_improvement_threshold": 0.1
    }
    arbiter = Arbiter(config)
    
    # 分数提升不足的情况
    history = [
        {"status": "failed", "score": 0.7}
    ]
    validation = ValidationResult(
        passed=False,
        score=0.75,  # 提升了0.05，小于阈值0.1
        issues=["测试错误"],
        suggestion="修改参数"
    )
    
    decision = arbiter.decide_next_action(history, validation)
    assert decision == "fallback"


def test_arbiter_decide_next_action_retry():
    """测试Arbiter decide_next_action方法 - 重试情况"""
    config = {
        "max_retries": 3,
        "score_improvement_threshold": 0.1
    }
    arbiter = Arbiter(config)
    
    # 其他情况，应该重试
    history = [
        {"status": "failed", "score": 0.6}
    ]
    validation = ValidationResult(
        passed=False,
        score=0.75,  # 提升了0.15，大于阈值0.1
        issues=["测试错误"],
        suggestion="修改参数"
    )
    
    decision = arbiter.decide_next_action(history, validation)
    assert decision == "retry"


def test_adaptive_iterator_initialization():
    """测试AdaptiveIterator初始化"""
    config = {
        "param_space": {
            "param1": {
                "type": "continuous",
                "min": 0.1,
                "max": 1.0,
                "default": 0.5
            },
            "param2": {
                "type": "discrete",
                "values": [1, 2, 3],
                "default": 2
            }
        }
    }
    iterator = AdaptiveIterator(config)
    assert iterator is not None
    assert "param1" in iterator.param_space
    assert "param2" in iterator.param_space


def test_adaptive_iterator_add_execution_record():
    """测试AdaptiveIterator add_execution_record方法"""
    config = {
        "param_space": {}
    }
    iterator = AdaptiveIterator(config)
    
    # 添加执行记录
    record = {
        "params": {"param1": 0.5},
        "validation": {
            "score": 0.8,
            "passed": True
        }
    }
    iterator.add_execution_record(record)
    
    history = iterator.get_execution_history()
    assert len(history) == 1
    assert history[0]["params"]["param1"] == 0.5


def test_adaptive_iterator_get_next_params():
    """测试AdaptiveIterator get_next_params方法"""
    config = {
        "param_space": {
            "param1": {
                "type": "continuous",
                "min": 0.1,
                "max": 1.0,
                "default": 0.5
            }
        }
    }
    iterator = AdaptiveIterator(config)
    
    # 获取下一组参数
    params = iterator.get_next_params()
    assert "param1" in params
    assert 0.1 <= params["param1"] <= 1.0


def test_adaptive_iterator_define_param_space():
    """测试AdaptiveIterator define_param_space方法"""
    config = {}
    iterator = AdaptiveIterator(config)
    
    # 定义参数空间
    param_space = {
        "param1": {
            "type": "continuous",
            "min": 0.1,
            "max": 1.0,
            "default": 0.5
        }
    }
    iterator.define_param_space(param_space)
    
    assert "param1" in iterator.param_space
    assert iterator.param_space["param1"]["default"] == 0.5


def test_adaptive_iterator_adjust_params():
    """测试AdaptiveIterator _adjust_params方法"""
    config = {
        "param_space": {
            "param1": {
                "type": "continuous",
                "min": 0.1,
                "max": 1.0,
                "default": 0.5
            }
        }
    }
    iterator = AdaptiveIterator(config)
    
    # 测试参数调整
    current_params = {"param1": 0.3}
    validation = {
        "score": 0.6,
        "passed": False
    }
    
    # 由于分数较低，应该向默认值靠拢
    new_params = iterator._adjust_params(current_params, validation)
    assert "param1" in new_params
    assert new_params["param1"] > 0.3  # 向默认值0.5靠拢
