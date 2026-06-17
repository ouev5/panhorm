#!/usr/bin/env python3

import argparse
import json
from src.rag_retriever import RAGRetriever


def main():
    """构建RAG索引"""
    parser = argparse.ArgumentParser(description='构建RAG索引')
    parser.add_argument('--documents', '-d', required=True, help='文档目录或文件路径')
    parser.add_argument('--output', '-o', default='./rag_index', help='索引输出目录')
    parser.add_argument('--report', '-r', default='rag_build_report.json', help='构建报告输出文件')
    
    args = parser.parse_args()
    
    print(f"Building RAG index from: {args.documents}")
    print(f"Output directory: {args.output}")
    
    # 初始化RAG检索器
    retriever = RAGRetriever(persist_directory=args.output)
    
    # 构建索引
    report = retriever.build_index(args.documents)
    
    # 输出构建报告
    print("\n=== Build Report ===")
    for key, value in report.items():
        print(f"{key}: {value}")
    
    # 保存报告到文件
    with open(args.report, 'w', encoding='utf-8') as f:
        json.dump(report, f, ensure_ascii=False, indent=2)
    
    print(f"\nReport saved to: {args.report}")
    print("Index build completed successfully!")


if __name__ == "__main__":
    main()
