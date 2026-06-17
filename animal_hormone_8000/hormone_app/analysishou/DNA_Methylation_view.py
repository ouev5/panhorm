import os
import json
import uuid
import logging
from django.shortcuts import render, redirect
from django.http import JsonResponse, HttpResponse, Http404
from django.views.decorators.csrf import csrf_exempt
from django.contrib.auth.decorators import login_required
from django.utils.decorators import method_decorator
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from django.conf import settings
from django.db import transaction
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import threading
import time
import tempfile
import zipfile
import csv
import io
from io import StringIO
import traceback
from pathlib import Path

# 导入你的模型
from ..models import DNAMethylationTask  # 根据实际路径调整

# 创建日志目录
LOG_DIR = Path('logs')
LOG_DIR.mkdir(exist_ok=True)

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(LOG_DIR / 'dna_methylation.log', encoding='utf-8'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger('dna_methylation')

# 创建详细的请求日志器
request_logger = logging.getLogger('dna_methylation.request')
request_logger.propagate = False
request_logger.setLevel(logging.INFO)
request_handler = logging.FileHandler(LOG_DIR / 'dna_methylation_requests.log', encoding='utf-8')
request_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
request_logger.addHandler(request_handler)

# 创建分析日志器
analysis_logger = logging.getLogger('dna_methylation.analysis')
analysis_logger.propagate = False
analysis_logger.setLevel(logging.INFO)
analysis_handler = logging.FileHandler(LOG_DIR / 'dna_methylation_analysis.log', encoding='utf-8')
analysis_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
analysis_logger.addHandler(analysis_handler)

@csrf_exempt
def test_dna_url(request, task_id=None):
    """测试DNA甲基化URL路由"""
    logger.info(f"测试URL访问 - 方法: {request.method}, 任务ID: {task_id}, 完整路径: {request.path}")
    return JsonResponse({
        'message': 'URL测试成功',
        'task_id': task_id,
        'path': request.path,
        'method': request.method
    })

class TaskLogger:
    """任务专用日志器"""
    
    def __init__(self, task_id):
        self.task_id = task_id
        self.task_logger = logging.getLogger(f'dna_methylation.task.{task_id}')
        self.task_logger.propagate = False
        self.task_logger.setLevel(logging.INFO)
        
        # 为每个任务创建单独的日志文件
        task_log_dir = LOG_DIR / 'tasks'
        task_log_dir.mkdir(exist_ok=True)
        
        task_handler = logging.FileHandler(
            task_log_dir / f'task_{task_id}.log', 
            encoding='utf-8'
        )
        task_handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s'
        ))
        self.task_logger.addHandler(task_handler)
    
    def info(self, message):
        self.task_logger.info(f"[{self.task_id}] {message}")
    
    def error(self, message):
        self.task_logger.error(f"[{self.task_id}] {message}")
    
    def warning(self, message):
        self.task_logger.warning(f"[{self.task_id}] {message}")
    
    def debug(self, message):
        self.task_logger.debug(f"[{self.task_id}] {message}")

class DNAMethylationAnalyzer:
    """DNA甲基化分析器"""
    
    def __init__(self, task_id):
        self.task_id = task_id
        self.logger = TaskLogger(task_id)
        self.status = "initialized"
        self.progress = 0
        self.results = {}
        self.logger.info("DNA甲基化分析器初始化")
    
    def process_methylation_data(self, file_path, params):
        """处理甲基化数据"""
        try:
            self.logger.info("开始处理甲基化数据")
            self.status = "processing"
            self.progress = 10
            
            # 记录参数
            self.logger.info(f"分析参数: {json.dumps(params, indent=2)}")
            
            # 读取数据文件
            self.logger.info(f"读取数据文件: {file_path}")
            if file_path.endswith('.csv') or file_path.endswith('.tsv'):
                delimiter = '\t' if file_path.endswith('.tsv') else ','
                df = pd.read_csv(file_path, sep=delimiter)
                self.logger.info(f"读取CSV/TSV文件成功，数据形状: {df.shape}")
                self.logger.info(f"列名: {df.columns.tolist()}")
            elif file_path.endswith('.bed'):
                # BED格式处理 - 动态处理列数
                try:
                    # 先尝试读取前几行来判断列数
                    with open(file_path, 'r') as f:
                        first_line = f.readline().strip()
                        col_count = len(first_line.split('\t'))
                    
                    self.logger.info(f"BED文件检测到 {col_count} 列")
                    
                    # 根据列数读取数据
                    if col_count >= 5:
                        # 读取所有列
                        df = pd.read_csv(file_path, sep='\t', header=None)
                        self.logger.info(f"读取BED文件成功，数据形状: {df.shape}")
                        
                        # 动态设置列名
                        if col_count == 5:
                            df.columns = ['chromosome', 'start', 'end', 'methylation_level', 'coverage']
                        elif col_count == 6:
                            df.columns = ['chromosome', 'start', 'end', 'name', 'methylation_level', 'coverage']
                        elif col_count == 7:
                            # 7列的BED格式：chr, start, end, name, score, strand, something
                            # 假设第5列是甲基化水平，第6列是覆盖率
                            df.columns = ['chromosome', 'start', 'end', 'name', 'methylation_level', 'coverage', 'strand']
                        elif col_count == 8:
                            df.columns = ['chromosome', 'start', 'end', 'name', 'methylation_level', 'coverage', 'strand', 'thickStart']
                        else:
                            # 对于更多列的情况，只取前7列
                            df = df.iloc[:, :7]
                            df.columns = ['chromosome', 'start', 'end', 'name', 'methylation_level', 'coverage', 'strand']
                        
                        self.logger.info(f"列名已设置为: {df.columns.tolist()}")
                        
                        # 检查必需的列是否存在
                        required_cols = ['chromosome', 'start', 'end', 'methylation_level', 'coverage']
                        missing_cols = [col for col in required_cols if col not in df.columns]
                        if missing_cols:
                            self.logger.warning(f"缺少列: {missing_cols}，将使用默认值")
                            # 为缺少的列添加默认值
                            for col in missing_cols:
                                if col == 'methylation_level':
                                    df[col] = 0.5  # 默认甲基化水平
                                elif col == 'coverage':
                                    df[col] = 10   # 默认覆盖率
                    else:
                        error_msg = f"BED文件列数不足: {col_count}列，至少需要5列"
                        self.logger.error(error_msg)
                        raise ValueError(error_msg)
                        
                except Exception as e:
                    error_msg = f"读取BED文件失败: {str(e)}"
                    self.logger.error(error_msg)
                    raise ValueError(error_msg)
            else:
                error_msg = f"不支持的文件格式: {file_path}"
                self.logger.error(error_msg)
                raise ValueError(error_msg)
            
            # 在质量控制之前，确保数值列是正确的类型
            self.logger.info("确保数值列类型正确")
            numeric_columns = ['start', 'end', 'methylation_level', 'coverage']
            for col in numeric_columns:
                if col in df.columns:
                    try:
                        # 记录转换前的类型
                        original_dtype = str(df[col].dtype)
                        sample_value = df[col].iloc[0] if len(df) > 0 else 'N/A'
                        self.logger.info(f"转换列 '{col}' - 原始类型: {original_dtype}, 示例值: {sample_value}")
                        
                        # 转换为数值类型，将错误转换为NaN
                        df[col] = pd.to_numeric(df[col], errors='coerce')
                        
                        # 检查转换后有多少NaN值
                        nan_count = df[col].isna().sum()
                        if nan_count > 0:
                            self.logger.warning(f"列 '{col}' 转换后有 {nan_count} 个NaN值 (共 {len(df)} 行)")
                            
                    except Exception as e:
                        self.logger.warning(f"转换列 '{col}' 时出错: {str(e)}")
                        # 如果转换失败，填充默认值
                        if col == 'methylation_level':
                            df[col] = 0.5
                        elif col == 'coverage':
                            df[col] = 10
            
            self.progress = 30
            self.logger.info(f"数据读取和类型转换完成，进度: {self.progress}%")
            
            # 质量控制
            self.logger.info("开始质量控制")
            filtered_df = self.quality_control(df, params)
            self.progress = 50
            self.logger.info(f"质量控制完成，过滤后数据形状: {filtered_df.shape}，进度: {self.progress}%")
            
            # 差异甲基化分析
            self.logger.info("开始差异甲基化分析")
            dmp_results = self.differential_methylation_analysis(filtered_df, params)
            self.results['dmp_results'] = dmp_results  # 供DMR区域合并复用，避免重复统计
            self.progress = 70
            self.logger.info(f"差异甲基化分析完成，找到 {len(dmp_results)} 个显著DMPs，进度: {self.progress}%")
            
            # 区域分析
            self.logger.info("开始区域分析")
            dmr_results = self.regional_analysis(filtered_df, params)
            self.progress = 85
            self.logger.info(f"区域分析完成，找到 {len(dmr_results)} 个DMRs，进度: {self.progress}%")
            
            # 功能富集分析
            self.logger.info("开始功能富集分析")
            enrichment_results = self.functional_enrichment_analysis(dmp_results, params)
            self.progress = 95
            self.logger.info(f"功能富集分析完成，进度: {self.progress}%")
            
            # 构建最终结果
            self.results = {
                'dmp_results': dmp_results,
                'dmr_results': dmr_results,
                'enrichment_results': enrichment_results,
                'summary_stats': self.calculate_summary_stats(filtered_df),
                'quality_metrics': self.calculate_quality_metrics(filtered_df)
            }
            
            self.progress = 100
            self.status = "completed"
            self.logger.info(f"分析完成，进度: {self.progress}%")
            
            # 记录摘要统计
            self.logger.info(f"分析摘要统计: {json.dumps(self.results['summary_stats'], indent=2)}")
            
            return self.results
        except Exception as e:
            error_msg = f"处理甲基化数据时出错: {str(e)}"
            self.logger.error(error_msg)
            self.logger.error(f"错误堆栈: {traceback.format_exc()}")
            self.status = "failed"
            self.progress = 0
            raise e
    
    def quality_control(self, df, params):
        """质量控制"""
        self.logger.info("执行质量控制")
        
        min_cpgs_per_sample = params.get('min_cpgs_per_sample', 50000)
        min_mean_coverage = params.get('min_mean_coverage', 5)
        max_meth_deviation = params.get('max_meth_deviation', 0.3)
        min_coverage_per_cpg = params.get('min_coverage_per_cpg', 3)
        min_samples_with_data = params.get('min_samples_with_data', 70)
        
        self.logger.info(f"质量控制参数: min_cpgs_per_sample={min_cpgs_per_sample}, "
                        f"min_coverage_per_cpg={min_coverage_per_cpg}")
        
        original_shape = df.shape
        self.logger.info(f"原始数据形状: {original_shape}")
        
        # 再次确保数值列是数值类型
        numeric_columns = ['coverage', 'methylation_level']
        for col in numeric_columns:
            if col in df.columns:
                if not pd.api.types.is_numeric_dtype(df[col]):
                    self.logger.warning(f"列 '{col}' 不是数值类型，正在转换...")
                    df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # 过滤低质量样本
        if 'sample' in df.columns:
            sample_counts = df.groupby('sample').size()
            self.logger.info(f"样本统计: 总样本数 {len(sample_counts)}")
            
            high_quality_samples = sample_counts[sample_counts >= min_cpgs_per_sample].index
            df = df[df['sample'].isin(high_quality_samples)]
            self.logger.info(f"样本过滤: 从 {len(sample_counts)} 个样本中保留 {len(high_quality_samples)} 个高质量样本")
        
        # 过滤低覆盖率CpG位点
        if 'coverage' in df.columns:
            # 确保coverage列是数值类型
            if not pd.api.types.is_numeric_dtype(df['coverage']):
                self.logger.warning("coverage列不是数值类型，强制转换...")
                df['coverage'] = pd.to_numeric(df['coverage'], errors='coerce')
            
            # 记录过滤前的统计信息
            before_count = len(df)
            if before_count > 0:
                self.logger.info(f"过滤前覆盖率统计: 最小值={df['coverage'].min()}, 最大值={df['coverage'].max()}, 平均值={df['coverage'].mean():.2f}")
            
            # 过滤低覆盖率
            df = df[df['coverage'].fillna(0) >= min_coverage_per_cpg]
            after_count = len(df)
            removed_count = before_count - after_count
            
            self.logger.info(f"覆盖率过滤: 移除了 {removed_count} 个低覆盖率CpG位点 (阈值: {min_coverage_per_cpg})")
            
            if after_count > 0:
                self.logger.info(f"过滤后覆盖率统计: 最小值={df['coverage'].min()}, 最大值={df['coverage'].max()}, 平均值={df['coverage'].mean():.2f}")
        
        # 过滤异常甲基化水平（通常应在0-1之间）
        if 'methylation_level' in df.columns:
            # 确保methylation_level列是数值类型
            if not pd.api.types.is_numeric_dtype(df['methylation_level']):
                self.logger.warning("methylation_level列不是数值类型，强制转换...")
                df['methylation_level'] = pd.to_numeric(df['methylation_level'], errors='coerce')
            
            before_count = len(df)
            # 甲基化水平应该在0-1之间
            df = df[(df['methylation_level'].fillna(-1) >= 0) & (df['methylation_level'].fillna(2) <= 1)]
            after_count = len(df)
            if after_count < before_count:
                self.logger.info(f"甲基化水平过滤: 移除了 {before_count - after_count} 个异常值 (范围: 0-1)")
            
            if after_count > 0:
                self.logger.info(f"甲基化水平统计: 最小值={df['methylation_level'].min():.3f}, 最大值={df['methylation_level'].max():.3f}, 平均值={df['methylation_level'].mean():.3f}")
        
        filtered_shape = df.shape
        self.logger.info(f"质量控制完成: 从 {original_shape} 过滤到 {filtered_shape}")
        self.logger.info(f"保留比例: {filtered_shape[0]/original_shape[0]*100:.1f}%")
        
        # 记录一些基本统计
        if not df.empty:
            self.logger.info(f"数据基本信息:")
            if 'chromosome' in df.columns:
                chrom_counts = df['chromosome'].value_counts()
                self.logger.info(f"  - 染色体分布: 共 {len(chrom_counts)} 条染色体")
                for chrom, count in chrom_counts.head(10).items():
                    self.logger.info(f"    {chrom}: {count} 个位点")
                if len(chrom_counts) > 10:
                    self.logger.info(f"    ... 和其他 {len(chrom_counts) - 10} 条染色体")
            
            if 'coverage' in df.columns:
                self.logger.info(f"  - 平均覆盖率: {df['coverage'].mean():.2f}")
                self.logger.info(f"  - 覆盖率中位数: {df['coverage'].median():.2f}")
            
            if 'methylation_level' in df.columns:
                self.logger.info(f"  - 平均甲基化水平: {df['methylation_level'].mean():.3f}")
                self.logger.info(f"  - 甲基化水平中位数: {df['methylation_level'].median():.3f}")
        
        return df
    
    def _bh_adjust(self, pvalues):
        """Benjamini-Hochberg FDR校正，返回与输入同序的q值。"""
        pvals = np.asarray([1.0 if pd.isna(p) else float(p) for p in pvalues], dtype=float)
        n = len(pvals)
        if n == 0:
            return []
        order = np.argsort(pvals)
        ranked = pvals[order]
        q = ranked * n / (np.arange(n) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        q = np.clip(q, 0, 1)
        out = np.empty(n, dtype=float)
        out[order] = q
        return out.tolist()

    def _load_metadata(self, metadata_path):
        """读取样本分组元数据，兼容CSV/TSV/TXT，标准化为 sample/group 两列。"""
        if not metadata_path or not os.path.exists(metadata_path):
            return None
        sep = '\t' if metadata_path.lower().endswith(('.tsv', '.txt')) else ','
        meta = pd.read_csv(metadata_path, sep=sep)
        lower = {str(c).strip().lower(): c for c in meta.columns}
        sample_col = next((lower[x] for x in ['sample', 'sample_id', 'sampleid', 'id'] if x in lower), meta.columns[0])
        group_col = next((lower[x] for x in ['group', 'condition', 'treatment', 'phenotype', 'class'] if x in lower), None)
        if not group_col or sample_col not in meta.columns:
            self.logger.warning('元数据缺少sample/group列，跳过分组差异分析')
            return None
        meta = meta[[sample_col, group_col]].dropna()
        meta.columns = ['sample', 'group']
        meta['sample'] = meta['sample'].astype(str)
        meta['group'] = meta['group'].astype(str)
        return meta

    def differential_methylation_analysis(self, df, params):
        """生产级DMP分析：基于真实样本分组执行逐CpG Welch t检验并做BH校正。"""
        self.logger.info("执行真实差异甲基化分析")
        pval_cutoff = float(params.get('pval_cutoff', params.get('detection_pval_thresh', 0.05)) or 0.05)
        beta_change_cutoff = float(params.get('beta_change_cutoff', 0.2) or 0.2)
        metadata_path = params.get('__metadata_file')
        max_sites = int(params.get('max_dmp_sites', 200000) or 200000)

        required = {'chromosome', 'start', 'methylation_level'}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"甲基化数据缺少必要列: {', '.join(sorted(missing))}")

        work = df.copy()
        if 'sample' not in work.columns:
            # 单样本数据不能做组间统计，返回按偏离全局均值排序的描述性结果（明确标注）
            self.logger.warning('数据无sample列，无法进行组间统计；输出描述性差异位点')
            baseline = float(work['methylation_level'].median()) if len(work) else 0.5
            work['beta_change'] = work['methylation_level'] - baseline
            work['p_value'] = 1.0
            work['adj_p_value'] = 1.0
            subset = work[work['beta_change'].abs() >= beta_change_cutoff].head(1000)
            return [{
                'chromosome': str(r.get('chromosome')),
                'position': int(r.get('start')),
                'gene': str(r.get('gene', r.get('gene_id', 'NA'))),
                'region': str(r.get('region', r.get('feature', 'NA'))),
                'beta_change': round(float(r.get('beta_change')), 6),
                'beta': round(float(r.get('beta_change')), 6),
                'p_value': 1.0,
                'adj_p_value': 1.0,
                'analysis_note': 'descriptive_only_no_sample_column'
            } for _, r in subset.iterrows()]

        work['sample'] = work['sample'].astype(str)
        meta = self._load_metadata(metadata_path)
        if meta is not None:
            work = work.merge(meta, on='sample', how='inner')
        elif 'group' not in work.columns:
            raise ValueError('差异甲基化分析需要元数据文件包含 sample/group 列，或甲基化数据自身包含 group 列')

        groups = [g for g in work['group'].dropna().unique().tolist()]
        if len(groups) < 2:
            raise ValueError('至少需要两个实验分组才能进行差异甲基化分析')
        ref = str(params.get('reference_group') or groups[0])
        test_groups = [g for g in groups if str(g) != ref]
        test = str(test_groups[0]) if test_groups else str(groups[1])
        self.logger.info(f"DMP比较: {test} vs {ref}; 阈值 beta={beta_change_cutoff}, FDR={pval_cutoff}")

        # 限制超大数据的站点数，避免Web进程被拖死；保留覆盖率最高/去重后的前N个CpG
        key_cols = ['chromosome', 'start']
        if len(work[key_cols].drop_duplicates()) > max_sites:
            self.logger.warning(f"CpG位点超过{max_sites}，按覆盖率优先截断以保障服务稳定")
            sort_col = 'coverage' if 'coverage' in work.columns else 'methylation_level'
            keep_keys = work.groupby(key_cols)[sort_col].mean().sort_values(ascending=False).head(max_sites).reset_index()[key_cols]
            work = work.merge(keep_keys, on=key_cols, how='inner')

        from scipy import stats
        records = []
        pvals = []
        grouped = work.groupby(key_cols, sort=False)
        for (chrom, pos), sub in grouped:
            a = pd.to_numeric(sub.loc[sub['group'].astype(str) == ref, 'methylation_level'], errors='coerce').dropna().values
            b = pd.to_numeric(sub.loc[sub['group'].astype(str) == test, 'methylation_level'], errors='coerce').dropna().values
            if len(a) < 2 or len(b) < 2:
                continue
            beta_change = float(np.mean(b) - np.mean(a))
            try:
                pval = float(stats.ttest_ind(b, a, equal_var=False, nan_policy='omit').pvalue)
            except Exception:
                pval = 1.0
            if not np.isfinite(pval):
                pval = 1.0
            row0 = sub.iloc[0]
            rec = {
                'chromosome': str(chrom),
                'position': int(pos),
                'gene': str(row0.get('gene', row0.get('gene_id', 'NA'))),
                'region': str(row0.get('region', row0.get('feature', 'NA'))),
                'beta_change': beta_change,
                'beta': beta_change,
                'p_value': pval,
                'mean_ref': float(np.mean(a)),
                'mean_test': float(np.mean(b)),
                'n_ref': int(len(a)),
                'n_test': int(len(b)),
                'comparison': f'{test}_vs_{ref}'
            }
            records.append(rec)
            pvals.append(pval)

        qvals = self._bh_adjust(pvals)
        for rec, q in zip(records, qvals):
            rec['adj_p_value'] = float(q)
        sig = [r for r in records if abs(r['beta_change']) >= beta_change_cutoff and r['adj_p_value'] <= pval_cutoff]
        sig.sort(key=lambda r: (r['adj_p_value'], -abs(r['beta_change'])))
        sig = sig[:5000]
        for r in sig:
            for k in ['beta_change', 'beta', 'p_value', 'adj_p_value', 'mean_ref', 'mean_test']:
                r[k] = round(float(r[k]), 8)
        self.logger.info(f"DMP真实统计完成: 测试位点 {len(records)}, 显著 {len(sig)}")
        return sig

    def regional_analysis(self, df, params):
        """生产级DMR分析：把显著DMP按染色体距离合并成差异甲基化区域。"""
        self.logger.info("执行真实区域DMR分析")
        max_gap = int(params.get('max_gap', 500) or 500)
        min_cpgs_dmr = int(params.get('min_cpgs_dmr', 3) or 3)
        min_region_length = int(params.get('min_region_length', 50) or 50)
        # 复用DMP真实统计结果，确保DMP/DMR一致
        dmp_results = self.results.get('dmp_results') if isinstance(self.results, dict) else None
        if not dmp_results:
            dmp_results = self.differential_methylation_analysis(df, params)
        if not dmp_results:
            self.logger.info("无显著DMP，DMR结果为空")
            return []
        dmp_df = pd.DataFrame(dmp_results)
        if dmp_df.empty:
            return []
        dmp_df['position'] = pd.to_numeric(dmp_df['position'], errors='coerce')
        dmp_df = dmp_df.dropna(subset=['chromosome', 'position']).sort_values(['chromosome', 'position'])
        regions = []
        for chrom, sub in dmp_df.groupby('chromosome', sort=False):
            current = []
            last_pos = None
            for _, row in sub.iterrows():
                pos = int(row['position'])
                if current and last_pos is not None and pos - last_pos > max_gap:
                    regions.append(current)
                    current = []
                current.append(row.to_dict())
                last_pos = pos
            if current:
                regions.append(current)
        dmrs = []
        for regs in regions:
            if len(regs) < min_cpgs_dmr:
                continue
            starts = [int(r['position']) for r in regs]
            start, end = min(starts), max(starts)
            if end - start + 1 < min_region_length:
                continue
            qvals = [float(r.get('adj_p_value', 1.0)) for r in regs]
            betas = [float(r.get('beta_change', r.get('beta', 0.0))) for r in regs]
            genes = [r.get('gene') for r in regs if r.get('gene') and r.get('gene') != 'NA']
            dmrs.append({
                'chromosome': str(regs[0]['chromosome']),
                'start': int(start),
                'end': int(end),
                'gene': ','.join(sorted(set(map(str, genes)))[:5]) if genes else 'NA',
                'cpg_count': int(len(regs)),
                'num_cpgs': int(len(regs)),
                'mean_beta_change': round(float(np.mean(betas)), 8),
                'beta': round(float(np.mean(betas)), 8),
                'p_value': round(float(min([float(r.get('p_value', 1.0)) for r in regs])), 8),
                'adj_p_value': round(float(min(qvals)), 8)
            })
        dmrs.sort(key=lambda x: (x['adj_p_value'], -x['cpg_count']))
        self.logger.info(f"DMR真实区域合并完成: {len(dmrs)} 个区域")
        return dmrs[:1000]

    def functional_enrichment_analysis(self, dmp_results, params):
        """轻量级真实富集：基于DMP关联基因做超几何检验；无本地注释库时返回可审计的基因摘要。"""
        self.logger.info("执行功能富集分析")
        enrichment_db = params.get('enrichment_db', ['go-bp', 'kegg']) or ['go-bp', 'kegg']
        genes = sorted({str(x.get('gene')) for x in (dmp_results or []) if x.get('gene') and str(x.get('gene')) != 'NA'})
        results = {}
        if not genes:
            return {db: [] for db in enrichment_db}
        # 内置常见生物过程关键词集；生产环境没有外部注释库时避免随机结果
        builtin_sets = {
            'hormone_response': {'ESR1','ESR2','AR','PGR','NR3C1','CYP19A1','STAR','LHCGR','FSHR','INSR'},
            'immune_response': {'IL6','TNF','IFNG','IL1B','CXCL8','CCL2','CD4','CD8A','NFKB1'},
            'cell_cycle': {'MKI67','PCNA','CDK1','CCND1','CCNE1','TOP2A','AURKA','BUB1'},
            'apoptosis': {'BAX','BCL2','CASP3','CASP8','TP53','FAS','BAD'},
            'stress_response': {'HSP90AA1','HSPA1A','HSPB1','SOD1','CAT','GPX1'}
        }
        from scipy.stats import hypergeom
        universe = max(20000, len(genes) * 20)
        for db in enrichment_db:
            terms = []
            for term, term_genes in builtin_sets.items():
                overlap = sorted(set(genes) & term_genes)
                if not overlap:
                    continue
                M, n, N, k = universe, len(term_genes), len(genes), len(overlap)
                pval = float(hypergeom.sf(k-1, M, n, N))
                terms.append({
                    'term_id': f'{str(db).upper()}:{term}',
                    'term_name': term.replace('_', ' ').title(),
                    'p_value': pval,
                    'adj_p_value': pval,
                    'gene_count': k,
                    'genes': overlap,
                    'gene_ratio': round(k / max(N, 1), 6),
                    'background_ratio': round(n / M, 6)
                })
            qvals = self._bh_adjust([t['p_value'] for t in terms])
            for t, q in zip(terms, qvals):
                t['adj_p_value'] = round(float(q), 8)
                t['p_value'] = round(float(t['p_value']), 8)
            terms.sort(key=lambda x: x['adj_p_value'])
            results[db] = terms[:50]
        self.logger.info(f"功能富集分析完成: 输入基因 {len(genes)}，富集项 {sum(len(v) for v in results.values())}")
        return results

    def calculate_summary_stats(self, df):
        """计算摘要统计"""
        self.logger.info("计算摘要统计")
        
        stats = {
            'total_cpgs': len(df),
            'samples': df['sample'].nunique() if 'sample' in df.columns else 1,
            'avg_methylation': float(df['methylation_level'].mean()) if 'methylation_level' in df.columns and len(df) > 0 else 0.5,
            'std_methylation': float(df['methylation_level'].std()) if 'methylation_level' in df.columns and len(df) > 1 else 0.1,
            'avg_coverage': float(df['coverage'].mean()) if 'coverage' in df.columns and len(df) > 0 else 10.0,
            'std_coverage': float(df['coverage'].std()) if 'coverage' in df.columns and len(df) > 1 else 2.0
        }
        
        self.logger.info(f"摘要统计: {json.dumps(stats, indent=2)}")
        return stats
    
    def calculate_quality_metrics(self, df):
        """计算质量指标"""
        self.logger.info("计算质量指标")
        
        metrics = {
            'coverage_range': (
                float(df['coverage'].min()) if 'coverage' in df.columns and len(df) > 0 else 0, 
                float(df['coverage'].max()) if 'coverage' in df.columns and len(df) > 0 else 100
            ),
            'methylation_range': (
                float(df['methylation_level'].min()) if 'methylation_level' in df.columns and len(df) > 0 else 0, 
                float(df['methylation_level'].max()) if 'methylation_level' in df.columns and len(df) > 0 else 1
            ),
            'missing_values': int(df.isnull().sum().sum()),
            'duplicate_positions': int(df.duplicated(subset=['chromosome', 'start', 'end']).sum()) if 'chromosome' in df.columns else 0
        }
        
        self.logger.info(f"质量指标: {json.dumps(metrics, indent=2)}")
        return metrics

def dna_methylation_analysis_page(request):
    """DNA甲基化分析页面"""
    logger.info(f"访问DNA甲基化分析页面 - IP: {request.META.get('REMOTE_ADDR')}")
    return render(request, 'dna_methylation_analysis.html')

@csrf_exempt
def dna_create_task(request):
    """创建DNA甲基化分析任务"""
    try:
        request_logger.info(f"创建任务请求 - 方法: {request.method}, IP: {request.META.get('REMOTE_ADDR')}")
        
        if request.method == 'POST':
            task_id = str(uuid.uuid4())
            
            # 创建数据库记录，确保files字段初始化为空列表
            task = DNAMethylationTask.objects.create(
                task_id=task_id,
                status='created',
                progress=0,
                params={},
                files=[],  # 初始化为空列表
                error=''
            )
            
            logger.info(f"创建DNA甲基化分析任务成功 - 任务ID: {task_id}")
            logger.info(f"任务创建详情: ID={task_id}, 状态={task.status}, 文件数={len(task.files) if task.files else 0}")
            request_logger.info(f"任务创建成功 - 任务ID: {task_id}")
            
            return JsonResponse({
                'task_id': task_id,
                'status': task.status,
                'message': '任务创建成功'
            })
        else:
            logger.warning(f"非法请求方法 - 方法: {request.method}")
            return JsonResponse({'error': 'Only POST method allowed'}, status=405)
            
    except Exception as e:
        error_msg = f"创建任务失败: {str(e)}"
        logger.error(error_msg)
        logger.error(f"错误堆栈: {traceback.format_exc()}")
        return JsonResponse({'error': error_msg}, status=500)

def detect_file_type(file_content, filename):
    """智能检测文件类型"""
    filename_lower = filename.lower()
    
    # 根据扩展名初步判断
    if filename_lower.endswith(('.bed', '.csv', '.tsv', '.txt')):
        # 检查文件内容是否为甲基化数据
        if 'chr' in file_content.lower() or any(keyword in file_content.lower() 
                for keyword in ['start', 'end', 'methyl', 'coverage', 'cg', 'position']):
            return 'methyl-data'
        else:
            # 检查是否是元数据
            if any(keyword in file_content.lower() 
                   for keyword in ['sample', 'group', 'treatment', 'condition', 'batch']):
                return 'metadata'
            else:
                # 默认为甲基化数据
                return 'methyl-data'
    elif filename_lower.endswith(('.gff', '.gtf')):
        return 'regions'
    elif filename_lower.endswith(('.gmt', '.gene', '.geneset')):
        return 'custom-genesets'
    else:
        return 'unknown'

@csrf_exempt
def dna_upload_file(request, task_id):
    """上传DNA甲基化分析文件"""
    try:
        request_logger.info(f"上传文件请求 - 任务ID: {task_id}, 方法: {request.method}")
        
        if request.method == 'POST':
            # 从数据库获取任务 - 需要重新获取最新状态
            try:
                task = DNAMethylationTask.objects.get(task_id=task_id)
            except DNAMethylationTask.DoesNotExist:
                error_msg = f"任务不存在: {task_id}"
                logger.error(error_msg)
                return JsonResponse({'error': error_msg}, status=404)
            
            # 获取上传的文件
            uploaded_file = request.FILES.get('file')
            if not uploaded_file:
                error_msg = "未上传文件"
                logger.warning(error_msg)
                return JsonResponse({'error': error_msg}, status=400)
            
            # 获取文件类型
            file_type = request.POST.get('file_type', 'unknown')
            
            # 如果类型未知，智能检测文件类型
            if file_type == 'unknown':
                # 读取文件前1KB内容进行检测
                file_content = uploaded_file.read(1024).decode('utf-8', errors='ignore')
                uploaded_file.seek(0)  # 重置文件指针
                file_type = detect_file_type(file_content, uploaded_file.name)
            
            logger.info(f"上传文件 - 任务ID: {task_id}, 文件名: {uploaded_file.name}, "
                       f"文件类型: {file_type}, 文件大小: {uploaded_file.size} bytes")
            
            # 保存文件
            file_path = default_storage.save(f'dna_methylation/{task_id}/{uploaded_file.name}', ContentFile(uploaded_file.read()))
            file_info = {
                'name': uploaded_file.name,
                'path': file_path,
                'type': file_type,
                'content_type': uploaded_file.content_type,
                'size': uploaded_file.size,
                'uploaded_at': datetime.now().isoformat()
            }
            
            # 使用原子操作更新任务记录 - 确保追加而不是覆盖
            with transaction.atomic():
                # 重新获取任务以确保我们有最新的数据
                task = DNAMethylationTask.objects.select_for_update().get(task_id=task_id)
                
                # 获取现有文件列表
                current_files = task.files if task.files else []
                
                # 检查是否已存在同名文件
                existing_files = [f for f in current_files if f.get('name') == uploaded_file.name]
                if existing_files:
                    # 更新现有文件信息
                    for i, file_item in enumerate(current_files):
                        if file_item.get('name') == uploaded_file.name:
                            current_files[i] = file_info
                            logger.info(f"更新现有文件: {uploaded_file.name}")
                            break
                else:
                    # 追加新文件
                    current_files.append(file_info)
                    logger.info(f"添加新文件: {uploaded_file.name}")
                
                # 更新任务记录
                task.files = current_files
                task.status = 'files_uploaded'
                task.save()
            
            # 重新获取任务以返回最新状态
            task.refresh_from_db()
            
            logger.info(f"文件上传成功 - 任务ID: {task_id}")
            logger.info(f"文件信息: {json.dumps(file_info, indent=2)}")
            logger.info(f"当前所有文件 ({len(task.files)} 个): {[f.get('name') for f in task.files]}")
            request_logger.info(f"文件上传成功 - 任务ID: {task_id}, 文件名: {uploaded_file.name}, 类型: {file_type}")
            
            return JsonResponse({
                'message': 'File uploaded successfully',
                'file_info': {
                    'name': uploaded_file.name,
                    'size': uploaded_file.size,
                    'type': file_type,
                    'content_type': uploaded_file.content_type
                },
                'total_files': len(task.files) if task.files else 0,
                'file_list': [{'name': f.get('name'), 'type': f.get('type')} for f in task.files]
            })
        else:
            logger.warning(f"非法请求方法 - 方法: {request.method}")
            return JsonResponse({'error': 'Only POST method allowed'}, status=405)
            
    except Exception as e:
        error_msg = f"上传文件失败: {str(e)}"
        logger.error(f"任务ID: {task_id}, {error_msg}")
        logger.error(f"错误堆栈: {traceback.format_exc()}")
        return JsonResponse({'error': error_msg}, status=500)

@csrf_exempt
def dna_run_analysis(request, task_id):
    """运行DNA甲基化分析"""
    try:
        request_logger.info(f"运行分析请求 - 任务ID: {task_id}, 方法: {request.method}")
        
        if request.method == 'POST':
            # 从数据库获取任务
            try:
                task = DNAMethylationTask.objects.get(task_id=task_id)
            except DNAMethylationTask.DoesNotExist:
                error_msg = f"任务不存在: {task_id}"
                logger.error(error_msg)
                return JsonResponse({'error': error_msg}, status=404)
            
            # 检查是否有文件
            if not task.files or len(task.files) == 0:
                error_msg = f"任务 {task_id} 未上传文件"
                logger.error(error_msg)
                return JsonResponse({'error': error_msg}, status=400)
            
            # 获取参数
            try:
                params = json.loads(request.body.decode('utf-8'))
                logger.info(f"分析参数 - 任务ID: {task_id}, 参数: {json.dumps(params, indent=2)}")
            except json.JSONDecodeError as e:
                error_msg = f"参数JSON解析失败: {str(e)}"
                logger.error(error_msg)
                return JsonResponse({'error': error_msg}, status=400)
            
            # 记录所有文件信息
            logger.info(f"分析文件信息 - 任务ID: {task_id}, 所有文件 ({len(task.files)} 个):")
            for i, file_info in enumerate(task.files):
                logger.info(f"  文件 {i+1}: {file_info.get('name')} (类型: {file_info.get('type')})")
            
            # 查找甲基化数据文件和元数据文件
            methyl_files = [f for f in task.files if f.get('type') == 'methyl-data']
            metadata_files = [f for f in task.files if f.get('type') == 'metadata']
            
            logger.info(f"甲基化数据文件: {[f.get('name') for f in methyl_files]}")
            logger.info(f"元数据文件: {[f.get('name') for f in metadata_files]}")
            
            # 如果没有明确标记为甲基化数据的文件，检查是否有可能是甲基化数据的文件
            if not methyl_files:
                potential_methyl_files = []
                for file_info in task.files:
                    file_name = file_info.get('name', '').lower()
                    file_type = file_info.get('type', 'unknown')
                    
                    # 如果文件类型是未知或者看起来像数据文件
                    if file_type == 'unknown' and file_name.endswith(('.bed', '.csv', '.tsv')):
                        potential_methyl_files.append(file_info)
                
                if potential_methyl_files:
                    logger.info(f"找到潜在的甲基化数据文件: {[f.get('name') for f in potential_methyl_files]}")
                    methyl_files = potential_methyl_files
            
            # 检查是否有甲基化数据文件
            if not methyl_files:
                error_msg = f"未找到甲基化数据文件。请上传BED、CSV或TSV格式的甲基化数据文件。现有文件: {[f.get('name') for f in task.files]}"
                logger.error(error_msg)
                
                # 更新任务状态为失败
                task.status = 'failed'
                task.error = error_msg
                task.save()
                
                return JsonResponse({'error': error_msg}, status=400)
            
            # 更新任务参数和状态
            task.params = params
            task.status = 'running'
            task.progress = 10
            task.error = ''
            task.save()
            
            logger.info(f"任务状态已更新为: {task.status}, 进度: {task.progress}%")
            
            # 在后台线程中运行分析
            def run_analysis_in_background():
                try:
                    task_logger = TaskLogger(task_id)
                    task_logger.info("开始后台分析任务")
                    
                    # 重新从数据库获取最新任务数据
                    try:
                        task = DNAMethylationTask.objects.get(task_id=task_id)
                    except DNAMethylationTask.DoesNotExist:
                        error_msg = f"任务不存在: {task_id}"
                        task_logger.error(error_msg)
                        raise ValueError(error_msg)
                    
                    # 选择甲基化数据文件（优先选择BED文件）
                    methyl_file = None
                    for file_info in methyl_files:
                        if file_info.get('name', '').lower().endswith('.bed'):
                            methyl_file = file_info
                            task_logger.info(f"选择BED文件: {methyl_file.get('name')}")
                            break
                    
                    # 如果没有BED文件，选择第一个文件
                    if not methyl_file and methyl_files:
                        methyl_file = methyl_files[0]
                        task_logger.info(f"选择第一个甲基化数据文件: {methyl_file.get('name')}")
                    
                    if not methyl_file:
                        error_msg = f"无法选择甲基化数据文件"
                        task_logger.error(error_msg)
                        raise ValueError(error_msg)
                    
                    full_file_path = os.path.join(settings.MEDIA_ROOT, methyl_file['path'])
                    # 生产级：把元数据文件路径传入分析器，用于真实分组差异甲基化统计
                    if metadata_files:
                        params = dict(params or {})
                        params['__metadata_file'] = os.path.join(settings.MEDIA_ROOT, metadata_files[0]['path'])
                        task_logger.info(f"使用元数据文件: {params['__metadata_file']}")
                    task_logger.info(f"使用甲基化数据文件: {full_file_path}")
                    
                    # 检查文件是否存在
                    if not os.path.exists(full_file_path):
                        error_msg = f"文件不存在: {full_file_path}"
                        task_logger.error(error_msg)
                        raise ValueError(error_msg)
                    
                    # 更新进度
                    task.progress = 20
                    task.save()
                    task_logger.info(f"分析进度更新: {task.progress}%")
                    
                    # 添加延迟，让前端有时间显示进度
                    time.sleep(1)
                    
                    # 运行分析
                    analyzer = DNAMethylationAnalyzer(task_id)
                    results = analyzer.process_methylation_data(full_file_path, params)
                    
                    # 更新任务结果
                    task.results = results
                    task.status = 'completed'
                    task.progress = 100
                    task.completed_at = datetime.now()
                    task.save()
                    
                    task_logger.info(f"分析完成 - 任务ID: {task_id}")
                    
                except Exception as e:
                    error_msg = f"后台分析失败: {str(e)}"
                    logger.error(f"任务ID: {task_id}, {error_msg}")
                    logger.error(f"错误堆栈: {traceback.format_exc()}")
                    
                    # 更新任务状态为失败
                    DNAMethylationTask.objects.filter(task_id=task_id).update(
                        status='failed',
                        error=str(e),
                        progress=0
                    )
            
            thread = threading.Thread(target=run_analysis_in_background, name=f"Analysis-{task_id}")
            thread.daemon = True
            thread.start()
            
            logger.info(f"分析任务已启动 - 任务ID: {task_id}, 线程: {thread.name}")
            request_logger.info(f"分析启动成功 - 任务ID: {task_id}")
            
            return JsonResponse({
                'message': 'Analysis started successfully', 
                'task_id': task_id,
                'status': 'running',
                'progress': task.progress,
                'methyl_files': [f.get('name') for f in methyl_files],
                'metadata_files': [f.get('name') for f in metadata_files]
            })
        else:
            logger.warning(f"非法请求方法 - 方法: {request.method}")
            return JsonResponse({'error': 'Only POST method allowed'}, status=405)
            
    except Exception as e:
        error_msg = f"启动分析失败: {str(e)}"
        logger.error(error_msg)
        logger.error(f"错误堆栈: {traceback.format_exc()}")
        return JsonResponse({'error': error_msg}, status=500)

@csrf_exempt
def dna_get_task_status(request, task_id):
    """获取DNA甲基化分析任务状态"""
    try:
        request_logger.debug(f"获取任务状态 - 任务ID: {task_id}")
        
        # 从数据库获取任务
        try:
            task = DNAMethylationTask.objects.get(task_id=task_id)
        except DNAMethylationTask.DoesNotExist:
            error_msg = f"任务不存在: {task_id}"
            logger.warning(error_msg)
            return JsonResponse({'error': error_msg}, status=404)
        
        status_info = {
            'task_id': task_id,
            'status': task.status,
            'progress': task.progress,
            'created_at': task.created_at.isoformat() if task.created_at else None,
            'completed_at': task.completed_at.isoformat() if task.completed_at else None,
            'files': task.files if task.files else [],
            'params': task.params if task.params else {},
            'error': task.error if task.error else None,
            'has_results': task.results is not None
        }
        
        # 如果是运行中状态，模拟进度更新（仅用于演示）
        if task.status == 'running' and task.progress < 95:
            current_progress = task.progress
            # 真实分析中应该由分析器更新进度
            # 这里只是演示，实际不应该在这里更新
            pass
        
        logger.debug(f"返回任务状态 - 任务ID: {task_id}, 状态: {task.status}, 进度: {task.progress}%")
        
        return JsonResponse(status_info)
        
    except Exception as e:
        error_msg = f"获取任务状态失败: {str(e)}"
        logger.error(error_msg)
        return JsonResponse({'error': error_msg}, status=500)

@csrf_exempt
def dna_get_task_results(request, task_id):
    """获取DNA甲基化分析结果"""
    try:
        request_logger.info(f"获取任务结果 - 任务ID: {task_id}")
        
        # 从数据库获取任务
        try:
            task = DNAMethylationTask.objects.get(task_id=task_id)
        except DNAMethylationTask.DoesNotExist:
            error_msg = f"任务不存在: {task_id}"
            logger.warning(error_msg)
            return JsonResponse({'error': error_msg}, status=404)
        
        if task.results is None:
            # 如果任务已完成但没有结果，返回特定状态
            if task.status == 'completed':
                error_msg = f"分析已完成但结果为空，请重新运行分析"
                logger.warning(error_msg)
                return JsonResponse({'error': error_msg}, status=404)
            else:
                status_msg = f"结果未就绪，当前状态: {task.status}"
                logger.info(status_msg)
                return JsonResponse({
                    'error': status_msg,
                    'status': task.status,
                    'progress': task.progress
                }, status=202)
        
        results = task.results
        
        # 转换字段名以匹配前端期望的格式
        processed_results = {
            'task_id': task_id,
            'total_cpgs': results.get('summary_stats', {}).get('total_cpgs', 0),
            'total_dmps': len(results.get('dmp_results', [])),
            'total_dmrs': len(results.get('dmr_results', [])),
            'associated_genes': len(set([dmp.get('gene', '') for dmp in results.get('dmp_results', []) if dmp.get('gene')])),
            'dmps': results.get('dmp_results', []),
            'dmrs': results.get('dmr_results', []),
            'summary_stats': results.get('summary_stats', {}),
            'quality_metrics': results.get('quality_metrics', {}),
            'enrichment_results': results.get('enrichment_results', {})
        }
        
        # 确保所有必要的字段都存在
        processed_results['total_cpgs'] = processed_results.get('total_cpgs', 0)
        processed_results['total_dmps'] = processed_results.get('total_dmps', 0)
        processed_results['total_dmrs'] = processed_results.get('total_dmrs', 0)
        processed_results['associated_genes'] = processed_results.get('associated_genes', 0)
        
        logger.info(f"返回任务结果 - 任务ID: {task_id}, "
                   f"总CpGs: {processed_results['total_cpgs']}, "
                   f"总DMPs: {processed_results['total_dmps']}, "
                   f"总DMRs: {processed_results['total_dmrs']}")
        
        return JsonResponse(processed_results)
        
    except Exception as e:
        error_msg = f"获取任务结果失败: {str(e)}"
        logger.error(error_msg)
        logger.error(f"错误堆栈: {traceback.format_exc()}")
        return JsonResponse({'error': error_msg}, status=500)

@csrf_exempt
def dna_download_result(request, result_id):
    """下载DNA甲基化分析结果"""
    try:
        task_id = result_id
        request_logger.info(f"下载结果 - 任务ID: {task_id}")
        
        # 从数据库获取任务
        try:
            task = DNAMethylationTask.objects.get(task_id=task_id)
        except DNAMethylationTask.DoesNotExist:
            error_msg = f"任务不存在: {task_id}"
            logger.warning(error_msg)
            return JsonResponse({'error': error_msg}, status=404)
        
        if task.results is None:
            error_msg = f"结果不可用: {task_id}"
            logger.warning(error_msg)
            return JsonResponse({'error': error_msg}, status=404)
        
        # 创建临时文件
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.json')
        try:
            with open(temp_file.name, 'w', encoding='utf-8') as f:
                json.dump(task.results, f, indent=2, ensure_ascii=False)
            
            with open(temp_file.name, 'rb') as f:
                response = HttpResponse(f.read(), content_type='application/json')
                response['Content-Disposition'] = f'attachment; filename=dna_methylation_results_{task_id}.json'
            
            logger.info(f"下载结果成功 - 任务ID: {task_id}, 文件大小: {os.path.getsize(temp_file.name)} bytes")
            return response
        finally:
            os.unlink(temp_file.name)
            
    except Exception as e:
        error_msg = f"下载结果失败: {str(e)}"
        logger.error(error_msg)
        logger.error(f"错误堆栈: {traceback.format_exc()}")
        return JsonResponse({'error': error_msg}, status=500)

@csrf_exempt
def dna_download_all_results(request, task_id):
    """下载所有DNA甲基化分析结果"""
    try:
        request_logger.info(f"下载所有结果 - 任务ID: {task_id}")
        
        # 从数据库获取任务
        try:
            task = DNAMethylationTask.objects.get(task_id=task_id)
        except DNAMethylationTask.DoesNotExist:
            error_msg = f"任务不存在: {task_id}"
            logger.warning(error_msg)
            return JsonResponse({'error': error_msg}, status=404)
        
        if task.results is None:
            error_msg = f"结果不可用: {task_id}"
            logger.warning(error_msg)
            return JsonResponse({'error': error_msg}, status=404)
        
        # 创建ZIP文件
        buffer = io.BytesIO()
        file_count = 0
        
        with zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            # 添加主要结果
            results_json = json.dumps(task.results, indent=2, ensure_ascii=False)
            zip_file.writestr(f'results_{task_id}.json', results_json)
            file_count += 1
            
            # 添加DMP结果CSV
            if 'dmp_results' in task.results:
                dmp_data = task.results['dmp_results']
                if dmp_data:
                    output = StringIO()
                    writer = csv.writer(output)
                    writer.writerow(['chromosome', 'position', 'gene', 'region', 'beta_change', 'p_value', 'adj_p_value'])
                    for dmp in dmp_data:
                        writer.writerow([
                            dmp['chromosome'], dmp['position'], dmp['gene'], dmp['region'],
                            dmp['beta_change'], dmp['p_value'], dmp['adj_p_value']
                        ])
                    zip_file.writestr(f'dmp_results_{task_id}.csv', output.getvalue())
                    file_count += 1
            
            # 添加DMR结果CSV
            if 'dmr_results' in task.results:
                dmr_data = task.results['dmr_results']
                if dmr_data:
                    output = StringIO()
                    writer = csv.writer(output)
                    writer.writerow(['chromosome', 'start', 'end', 'gene', 'cpg_count', 'mean_beta_change', 'p_value', 'adj_p_value'])
                    for dmr in dmr_data:
                        writer.writerow([
                            dmr['chromosome'], dmr['start'], dmr['end'], dmr['gene'],
                            dmr['cpg_count'], dmr['mean_beta_change'], dmr['p_value'], dmr['adj_p_value']
                        ])
                    zip_file.writestr(f'dmr_results_{task_id}.csv', output.getvalue())
                    file_count += 1
            
            # 添加摘要统计
            if 'summary_stats' in task.results:
                summary_json = json.dumps(task.results['summary_stats'], indent=2, ensure_ascii=False)
                zip_file.writestr(f'summary_stats_{task_id}.json', summary_json)
                file_count += 1
            
            # 添加质量指标
            if 'quality_metrics' in task.results:
                quality_json = json.dumps(task.results['quality_metrics'], indent=2, ensure_ascii=False)
                zip_file.writestr(f'quality_metrics_{task_id}.json', quality_json)
                file_count += 1
            
            # 添加富集分析结果
            if 'enrichment_results' in task.results:
                enrichment_json = json.dumps(task.results['enrichment_results'], indent=2, ensure_ascii=False)
                zip_file.writestr(f'enrichment_results_{task_id}.json', enrichment_json)
                file_count += 1
        
        buffer.seek(0)
        response = HttpResponse(buffer.getvalue(), content_type='application/zip')
        response['Content-Disposition'] = f'attachment; filename=dna_methylation_all_results_{task_id}.zip'
        
        logger.info(f"下载所有结果成功 - 任务ID: {task_id}, 文件数: {file_count}, ZIP大小: {buffer.tell()} bytes")
        request_logger.info(f"下载所有结果成功 - 任务ID: {task_id}")
        
        return response
        
    except Exception as e:
        error_msg = f"下载所有结果失败: {str(e)}"
        logger.error(error_msg)
        logger.error(f"错误堆栈: {traceback.format_exc()}")
        return JsonResponse({'error': error_msg}, status=500)

@csrf_exempt
@csrf_exempt
def dna_get_visualization_plot(request, task_id, plot_type):
    """获取真实结果驱动的甲基化可视化Plotly数据。"""
    try:
        request_logger.info(f"获取可视化图表 - 任务ID: {task_id}, 图表类型: {plot_type}")
        try:
            task = DNAMethylationTask.objects.get(task_id=task_id)
        except DNAMethylationTask.DoesNotExist:
            return JsonResponse({'error': f'任务不存在: {task_id}'}, status=404)
        if task.results is None:
            return JsonResponse({'error': f'结果不可用: {task_id}'}, status=404)
        results = task.results or {}
        dmps = results.get('dmp_results', []) or []
        dmrs = results.get('dmr_results', []) or []
        summary = results.get('summary_stats', {}) or {}
        quality = results.get('quality_metrics', {}) or {}

        if plot_type == 'methylation-distribution':
            avg = float(summary.get('avg_methylation', 0.5) or 0.5)
            std = max(float(summary.get('std_methylation', 0.1) or 0.1), 0.02)
            x = np.linspace(0, 1, 50)
            y = np.exp(-0.5*((x-avg)/std)**2)
            y = (y / y.sum() * int(summary.get('total_cpgs', 100) or 100)).astype(int).tolist()
            data = {'data':[{'x':[round(float(v),3) for v in x], 'y':y, 'type':'bar', 'name':'Estimated distribution'}], 'layout':{'title':'Methylation Beta Value Distribution','xaxis':{'title':'Beta value'},'yaxis':{'title':'CpG count'}}}
        elif plot_type == 'dmp-volcano':
            x = [float(d.get('beta_change', d.get('beta', 0))) for d in dmps]
            pvals = [max(float(d.get('adj_p_value', d.get('p_value', 1)) or 1), 1e-300) for d in dmps]
            y = [-float(np.log10(p)) for p in pvals]
            genes = [str(d.get('gene','')) for d in dmps]
            colors = ['#e11d48' if abs(a) >= 0.2 and b >= -np.log10(0.05) else '#64748b' for a,b in zip(x,y)]
            data = {'data':[{'x':x,'y':y,'type':'scatter','mode':'markers','marker':{'color':colors,'size':7},'text':genes,'name':'DMPs'}], 'layout':{'title':'Volcano Plot of Differentially Methylated Positions','xaxis':{'title':'Delta beta'},'yaxis':{'title':'-log10(FDR)'}}}
        elif plot_type in ['dmr-browser','regional-methylation','chromosomal-distribution']:
            if plot_type == 'chromosomal-distribution':
                counts = pd.Series([d.get('chromosome') for d in dmps]).value_counts().sort_index()
                data = {'data':[{'x':counts.index.astype(str).tolist(),'y':counts.astype(int).tolist(),'type':'bar','name':'DMP count'}], 'layout':{'title':'Chromosomal Distribution of DMPs','xaxis':{'title':'Chromosome'},'yaxis':{'title':'Count'}}}
            else:
                data = {'data':[{'x':[f"{d.get('chromosome')}:{d.get('start')}-{d.get('end')}" for d in dmrs[:50]], 'y':[float(d.get('mean_beta_change', d.get('beta',0))) for d in dmrs[:50]], 'type':'bar', 'name':'DMR delta beta'}], 'layout':{'title':'Differentially Methylated Regions','xaxis':{'title':'Region'},'yaxis':{'title':'Mean delta beta'}}}
        elif plot_type == 'pca-analysis':
            # 无原始矩阵时不伪造PCA，返回质量指标散点（真实统计）
            data = {'data':[{'x':[summary.get('avg_coverage',0)], 'y':[summary.get('avg_methylation',0)], 'type':'scatter','mode':'markers','marker':{'size':14},'name':'Dataset'}], 'layout':{'title':'Dataset QC Summary (coverage vs methylation)','xaxis':{'title':'Average coverage'},'yaxis':{'title':'Average methylation'}}}
        elif plot_type == 'enrichment-barplot':
            enrichment = results.get('enrichment_results', {}) or {}
            terms = []
            for arr in enrichment.values(): terms.extend(arr or [])
            terms = sorted(terms, key=lambda t: float(t.get('adj_p_value',1)))[:15]
            data = {'data':[{'x':[-np.log10(max(float(t.get('adj_p_value',1)),1e-300)) for t in terms], 'y':[t.get('term_name', t.get('term_id','')) for t in terms], 'type':'bar','orientation':'h','name':'Enrichment'}], 'layout':{'title':'Top Enriched Terms','xaxis':{'title':'-log10(FDR)'},'height':500}}
        elif plot_type == 'methylation-heatmap':
            top = dmps[:30]
            data = {'data':[{'z':[[float(d.get('mean_ref',0)), float(d.get('mean_test',0))] for d in top], 'x':['Reference','Test'], 'y':[f"{d.get('gene','NA')}@{d.get('chromosome')}:{d.get('position')}" for d in top], 'type':'heatmap','colorscale':'Viridis'}], 'layout':{'title':'Top DMP Mean Methylation Heatmap','height':700}}
        elif plot_type == 'meth-expression':
            data = {'data':[{'x':[float(d.get('beta_change',0)) for d in dmps[:200]], 'y':[-np.log10(max(float(d.get('adj_p_value',1)),1e-300)) for d in dmps[:200]], 'type':'scatter','mode':'markers','text':[d.get('gene','') for d in dmps[:200]]}], 'layout':{'title':'Methylation Effect vs Significance','xaxis':{'title':'Delta beta'},'yaxis':{'title':'-log10(FDR)'}}}
        else:
            return JsonResponse({'error': f'未知或未生成的图表类型: {plot_type}'}, status=404)
        return JsonResponse(data)
    except Exception as e:
        logger.error(f"获取图表数据失败: {str(e)}")
        logger.error(f"错误堆栈: {traceback.format_exc()}")
        return JsonResponse({'error': f'获取图表数据失败: {str(e)}'}, status=500)

@csrf_exempt
def dna_download_result_by_type(request, task_id):
    """按类型下载DNA甲基化分析结果"""
    try:
        request_logger.info(f"按类型下载结果 - 任务ID: {task_id}")
        
        # 从数据库获取任务
        try:
            task = DNAMethylationTask.objects.get(task_id=task_id)
        except DNAMethylationTask.DoesNotExist:
            error_msg = f"任务不存在: {task_id}"
            logger.warning(error_msg)
            return JsonResponse({'error': error_msg}, status=404)
        
        if task.results is None:
            error_msg = f"结果不可用: {task_id}"
            logger.warning(error_msg)
            return JsonResponse({'error': error_msg}, status=404)
        
        # 获取result_type参数
        result_type = 'all'
        if request.method == 'POST':
            try:
                request_data = json.loads(request.body.decode('utf-8'))
                result_type = request_data.get('result_type', 'all')
            except:
                result_type = 'all'
        else:
            result_type = request.GET.get('result_type', 'all')
        
        logger.info(f"下载类型: {result_type} - 任务ID: {task_id}")
        
        if result_type == 'dmps':
            # 下载DMP结果
            if 'dmp_results' not in task.results:
                error_msg = f"DMP结果不可用: {task_id}"
                logger.warning(error_msg)
                return JsonResponse({'error': error_msg}, status=404)
            
            output = StringIO()
            writer = csv.writer(output)
            writer.writerow(['chromosome', 'position', 'gene', 'region', 'beta_change', 'p_value', 'adj_p_value'])
            for dmp in task.results['dmp_results']:
                writer.writerow([
                    dmp['chromosome'], dmp['position'], dmp['gene'], dmp['region'],
                    dmp['beta_change'], dmp['p_value'], dmp['adj_p_value']
                ])
            
            response = HttpResponse(output.getvalue(), content_type='text/csv')
            response['Content-Disposition'] = f'attachment; filename=dna_methylation_dmp_results_{task_id}.csv'
            
            logger.info(f"DMP结果下载成功 - 任务ID: {task_id}, 记录数: {len(task.results['dmp_results'])}")
            return response
        
        elif result_type == 'dmrs':
            # 下载DMR结果
            if 'dmr_results' not in task.results:
                error_msg = f"DMR结果不可用: {task_id}"
                logger.warning(error_msg)
                return JsonResponse({'error': error_msg}, status=404)
            
            output = StringIO()
            writer = csv.writer(output)
            writer.writerow(['chromosome', 'start', 'end', 'gene', 'cpg_count', 'mean_beta_change', 'p_value', 'adj_p_value'])
            for dmr in task.results['dmr_results']:
                writer.writerow([
                    dmr['chromosome'], dmr['start'], dmr['end'], dmr['gene'],
                    dmr['cpg_count'], dmr['mean_beta_change'], dmr['p_value'], dmr['adj_p_value']
                ])
            
            response = HttpResponse(output.getvalue(), content_type='text/csv')
            response['Content-Disposition'] = f'attachment; filename=dna_methylation_dmr_results_{task_id}.csv'
            
            logger.info(f"DMR结果下载成功 - 任务ID: {task_id}, 记录数: {len(task.results['dmr_results'])}")
            return response
        
        elif result_type == 'enrichment':
            # 下载富集分析结果
            if 'enrichment_results' not in task.results:
                error_msg = f"富集结果不可用: {task_id}"
                logger.warning(error_msg)
                return JsonResponse({'error': error_msg}, status=404)
            
            output = StringIO()
            writer = csv.writer(output)
            writer.writerow(['term_id', 'term_name', 'p_value', 'adj_p_value', 'gene_count', 'gene_ratio', 'background_ratio'])
            
            total_terms = 0
            for db, terms in task.results['enrichment_results'].items():
                for term in terms:
                    writer.writerow([
                        term['term_id'], term['term_name'], term['p_value'], 
                        term['adj_p_value'], term['gene_count'], term['gene_ratio'], term['background_ratio']
                    ])
                    total_terms += 1
            
            response = HttpResponse(output.getvalue(), content_type='text/csv')
            response['Content-Disposition'] = f'attachment; filename=dna_methylation_enrichment_results_{task_id}.csv'
            
            logger.info(f"富集结果下载成功 - 任务ID: {task_id}, 数据库数: {len(task.results['enrichment_results'])}, 总项数: {total_terms}")
            return response
        
        else:
            # 默认下载所有结果
            logger.info(f"下载所有结果 - 任务ID: {task_id}")
            return dna_download_all_results(request, task_id)
            
    except Exception as e:
        error_msg = f"按类型下载结果失败: {str(e)}"
        logger.error(f"任务ID: {task_id}, 下载类型: {result_type}, {error_msg}")
        logger.error(f"错误堆栈: {traceback.format_exc()}")
        return JsonResponse({'error': error_msg}, status=500)

def cleanup_old_tasks():
    """清理旧任务（可选，可用于定时任务）"""
    try:
        logger.info("开始清理旧任务")
        current_time = datetime.now()
        
        # 清理24小时前的任务
        cutoff_time = current_time - timedelta(days=1)
        
        # 删除旧任务
        old_tasks = DNAMethylationTask.objects.filter(created_at__lt=cutoff_time)
        count = old_tasks.count()
        old_tasks.delete()
        
        logger.info(f"清理完成，移除 {count} 个旧任务")
        
    except Exception as e:
        logger.error(f"清理旧任务失败: {str(e)}")