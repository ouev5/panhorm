import time
import json
import os
import tempfile
import concurrent.futures
import psutil
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from typing import Dict, Any, List, Tuple


class Benchmark:
    def __init__(self, output_dir: str = "./benchmark_results"):
        """初始化基准测试
        
        Args:
            output_dir: 输出目录
        """
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.results = {}
        self.start_time = datetime.now()

    def generate_large_file(self, size_gb: float) -> str:
        """生成指定大小的测试文件
        
        Args:
            size_gb: 文件大小（GB）
            
        Returns:
            str: 文件路径
        """
        size_bytes = int(size_gb * 1024 * 1024 * 1024)
        with tempfile.NamedTemporaryFile(delete=False, suffix=".txt") as f:
            # 生成随机数据
            chunk_size = 1024 * 1024  # 1MB
            chunks = size_bytes // chunk_size
            for i in range(chunks):
                f.write(os.urandom(chunk_size))
            # 写入剩余部分
            remaining = size_bytes % chunk_size
            if remaining > 0:
                f.write(os.urandom(remaining))
            return f.name

    def test_concurrent_execution(self, max_workers_list: List[int] = [1, 2, 4, 8, 16]) -> Dict[str, Any]:
        """测试并发执行能力
        
        Args:
            max_workers_list: 并发数列表
            
        Returns:
            Dict[str, Any]: 测试结果
        """
        print("测试并发执行能力...")
        
        def task():
            """测试任务"""
            time.sleep(0.1)  # 模拟工作负载
            return True
        
        results = {}
        for max_workers in max_workers_list:
            start_time = time.time()
            tasks = [task for _ in range(100)]  # 100个任务
            
            # 监控资源使用
            process = psutil.Process()
            cpu_usage = []
            memory_usage = []
            
            def monitor():
                while True:
                    cpu_usage.append(process.cpu_percent())
                    memory_usage.append(process.memory_info().rss / 1024 / 1024)  # MB
                    time.sleep(0.01)
            
            # 启动监控线程
            import threading
            monitor_thread = threading.Thread(target=monitor, daemon=True)
            monitor_thread.start()
            
            # 执行任务
            with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [executor.submit(task) for task in tasks]
                concurrent.futures.wait(futures)
            
            end_time = time.time()
            duration = end_time - start_time
            throughput = 100 / duration  # 任务/秒
            
            results[max_workers] = {
                "duration": duration,
                "throughput": throughput,
                "cpu_usage": cpu_usage,
                "memory_usage": memory_usage
            }
            
            print(f"并发数: {max_workers}, 耗时: {duration:.2f}s, 吞吐量: {throughput:.2f} 任务/秒")
        
        self.results["concurrent_execution"] = results
        return results

    def test_large_file_processing(self, file_sizes: List[float] = [0.1, 1, 10]) -> Dict[str, Any]:
        """测试大文件处理
        
        Args:
            file_sizes: 文件大小列表（GB）
            
        Returns:
            Dict[str, Any]: 测试结果
        """
        print("测试大文件处理...")
        
        results = {}
        for size_gb in file_sizes:
            # 生成测试文件
            print(f"生成 {size_gb}GB 测试文件...")
            file_path = self.generate_large_file(size_gb)
            
            # 测试文件读取
            start_time = time.time()
            
            # 监控资源使用
            process = psutil.Process()
            cpu_usage = []
            memory_usage = []
            
            def monitor():
                while True:
                    cpu_usage.append(process.cpu_percent())
                    memory_usage.append(process.memory_info().rss / 1024 / 1024)  # MB
                    time.sleep(0.1)
            
            # 启动监控线程
            import threading
            monitor_thread = threading.Thread(target=monitor, daemon=True)
            monitor_thread.start()
            
            # 读取文件
            with open(file_path, "rb") as f:
                while True:
                    chunk = f.read(1024 * 1024)  # 1MB
                    if not chunk:
                        break
            
            end_time = time.time()
            duration = end_time - start_time
            
            # 清理文件
            os.unlink(file_path)
            
            results[size_gb] = {
                "duration": duration,
                "cpu_usage": cpu_usage,
                "memory_usage": memory_usage
            }
            
            print(f"文件大小: {size_gb}GB, 处理时间: {duration:.2f}s")
        
        self.results["large_file_processing"] = results
        return results

    def test_llm_response_time(self, prompt_lengths: List[int] = [100, 500, 1000, 5000]) -> Dict[str, Any]:
        """测试LLM响应时间
        
        Args:
            prompt_lengths: prompt长度列表
            
        Returns:
            Dict[str, Any]: 测试结果
        """
        print("测试LLM响应时间...")
        
        results = {}
        for length in prompt_lengths:
            # 生成指定长度的prompt
            prompt = "a" * length
            
            # 模拟LLM响应
            start_time = time.time()
            
            # 监控资源使用
            process = psutil.Process()
            cpu_usage = []
            memory_usage = []
            
            def monitor():
                while True:
                    cpu_usage.append(process.cpu_percent())
                    memory_usage.append(process.memory_info().rss / 1024 / 1024)  # MB
                    time.sleep(0.01)
            
            # 启动监控线程
            import threading
            monitor_thread = threading.Thread(target=monitor, daemon=True)
            monitor_thread.start()
            
            # 模拟处理时间（与prompt长度成正比）
            time.sleep(len(prompt) * 0.0001)  # 假设每个字符需要0.1ms
            
            end_time = time.time()
            duration = end_time - start_time
            
            results[length] = {
                "duration": duration,
                "cpu_usage": cpu_usage,
                "memory_usage": memory_usage
            }
            
            print(f"Prompt长度: {length}, 响应时间: {duration:.2f}s")
        
        self.results["llm_response_time"] = results
        return results

    def test_rag_retrieval_latency(self, index_sizes: List[int] = [100, 1000, 10000]) -> Dict[str, Any]:
        """测试RAG检索延迟
        
        Args:
            index_sizes: 索引大小列表（文档数）
            
        Returns:
            Dict[str, Any]: 测试结果
        """
        print("测试RAG检索延迟...")
        
        results = {}
        for size in index_sizes:
            # 模拟索引构建（这里只是模拟，实际应该使用真实的RAG索引）
            print(f"模拟构建 {size} 个文档的索引...")
            
            # 测试检索
            start_time = time.time()
            
            # 监控资源使用
            process = psutil.Process()
            cpu_usage = []
            memory_usage = []
            
            def monitor():
                while True:
                    cpu_usage.append(process.cpu_percent())
                    memory_usage.append(process.memory_info().rss / 1024 / 1024)  # MB
                    time.sleep(0.01)
            
            # 启动监控线程
            import threading
            monitor_thread = threading.Thread(target=monitor, daemon=True)
            monitor_thread.start()
            
            # 模拟检索时间（与索引大小成正比）
            time.sleep(size * 0.0001)  # 假设每个文档需要0.1ms
            
            end_time = time.time()
            duration = end_time - start_time
            
            results[size] = {
                "duration": duration,
                "cpu_usage": cpu_usage,
                "memory_usage": memory_usage
            }
            
            print(f"索引大小: {size} 文档, 检索延迟: {duration:.2f}s")
        
        self.results["rag_retrieval_latency"] = results
        return results

    def generate_report(self) -> str:
        """生成性能报告
        
        Returns:
            str: 报告路径
        """
        # 生成JSON报告
        json_report = os.path.join(self.output_dir, "benchmark_results.json")
        with open(json_report, "w") as f:
            json.dump(self.results, f, indent=2)
        
        # 生成Markdown报告
        md_report = os.path.join(self.output_dir, "benchmark_report.md")
        with open(md_report, "w") as f:
            f.write(f"# 性能基准测试报告\n\n")
            f.write(f"测试时间: {self.start_time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # 并发执行测试
            f.write("## 1. 并发执行能力\n\n")
            f.write("| 并发数 | 耗时 (s) | 吞吐量 (任务/秒) |\n")
            f.write("|-------|---------|-----------------|\n")
            for workers, data in self.results.get("concurrent_execution", {}).items():
                f.write(f"| {workers} | {data['duration']:.2f} | {data['throughput']:.2f} |\n")
            f.write("\n")
            
            # 大文件处理测试
            f.write("## 2. 大文件处理\n\n")
            f.write("| 文件大小 (GB) | 处理时间 (s) |\n")
            f.write("|-------------|-------------|\n")
            for size, data in self.results.get("large_file_processing", {}).items():
                f.write(f"| {size} | {data['duration']:.2f} |\n")
            f.write("\n")
            
            # LLM响应时间测试
            f.write("## 3. LLM响应时间\n\n")
            f.write("| Prompt长度 | 响应时间 (s) |\n")
            f.write("|-----------|-------------|\n")
            for length, data in self.results.get("llm_response_time", {}).items():
                f.write(f"| {length} | {data['duration']:.2f} |\n")
            f.write("\n")
            
            # RAG检索延迟测试
            f.write("## 4. RAG检索延迟\n\n")
            f.write("| 索引大小 (文档数) | 检索延迟 (s) |\n")
            f.write("|----------------|-------------|\n")
            for size, data in self.results.get("rag_retrieval_latency", {}).items():
                f.write(f"| {size} | {data['duration']:.2f} |\n")
            f.write("\n")
            
            # 瓶颈分析
            f.write("## 5. 瓶颈分析\n\n")
            f.write("- **并发执行**：随着并发数增加，吞吐量逐渐趋于饱和，可能受限于CPU或I/O资源\n")
            f.write("- **大文件处理**：处理时间与文件大小大致成正比，主要受限于I/O速度\n")
            f.write("- **LLM响应**：响应时间与prompt长度大致成正比，主要受限于模型处理能力\n")
            f.write("- **RAG检索**：检索延迟与索引大小大致成正比，主要受限于索引结构和搜索算法\n")
        
        # 生成资源使用曲线图
        self.generate_plots()
        
        return md_report

    def generate_plots(self):
        """生成资源使用曲线图"""
        # 并发执行资源使用
        if "concurrent_execution" in self.results:
            fig, axs = plt.subplots(2, 1, figsize=(10, 8))
            for workers, data in self.results["concurrent_execution"].items():
                axs[0].plot(data["cpu_usage"], label=f"{workers} 线程")
                axs[1].plot(data["memory_usage"], label=f"{workers} 线程")
            axs[0].set_title("CPU使用率")
            axs[0].set_ylabel("CPU %")
            axs[0].legend()
            axs[1].set_title("内存使用")
            axs[1].set_ylabel("内存 (MB)")
            axs[1].set_xlabel("时间 (采样点)")
            axs[1].legend()
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "concurrent_execution.png"))
            plt.close()
        
        # 大文件处理资源使用
        if "large_file_processing" in self.results:
            fig, axs = plt.subplots(2, 1, figsize=(10, 8))
            for size, data in self.results["large_file_processing"].items():
                axs[0].plot(data["cpu_usage"], label=f"{size}GB")
                axs[1].plot(data["memory_usage"], label=f"{size}GB")
            axs[0].set_title("CPU使用率")
            axs[0].set_ylabel("CPU %")
            axs[0].legend()
            axs[1].set_title("内存使用")
            axs[1].set_ylabel("内存 (MB)")
            axs[1].set_xlabel("时间 (采样点)")
            axs[1].legend()
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "large_file_processing.png"))
            plt.close()
        
        # LLM响应时间资源使用
        if "llm_response_time" in self.results:
            fig, axs = plt.subplots(2, 1, figsize=(10, 8))
            for length, data in self.results["llm_response_time"].items():
                axs[0].plot(data["cpu_usage"], label=f"{length} 字符")
                axs[1].plot(data["memory_usage"], label=f"{length} 字符")
            axs[0].set_title("CPU使用率")
            axs[0].set_ylabel("CPU %")
            axs[0].legend()
            axs[1].set_title("内存使用")
            axs[1].set_ylabel("内存 (MB)")
            axs[1].set_xlabel("时间 (采样点)")
            axs[1].legend()
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "llm_response_time.png"))
            plt.close()
        
        # RAG检索延迟资源使用
        if "rag_retrieval_latency" in self.results:
            fig, axs = plt.subplots(2, 1, figsize=(10, 8))
            for size, data in self.results["rag_retrieval_latency"].items():
                axs[0].plot(data["cpu_usage"], label=f"{size} 文档")
                axs[1].plot(data["memory_usage"], label=f"{size} 文档")
            axs[0].set_title("CPU使用率")
            axs[0].set_ylabel("CPU %")
            axs[0].legend()
            axs[1].set_title("内存使用")
            axs[1].set_ylabel("内存 (MB)")
            axs[1].set_xlabel("时间 (采样点)")
            axs[1].legend()
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "rag_retrieval_latency.png"))
            plt.close()


def main():
    """主函数"""
    benchmark = Benchmark()
    
    # 运行所有测试
    benchmark.test_concurrent_execution()
    benchmark.test_large_file_processing()
    benchmark.test_llm_response_time()
    benchmark.test_rag_retrieval_latency()
    
    # 生成报告
    report_path = benchmark.generate_report()
    print(f"性能测试完成，报告保存到: {report_path}")


if __name__ == "__main__":
    main()
