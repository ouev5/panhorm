# /www/wwwroot/default/animal_hormone/hormone_app/management/commands/predict_by_network_fixed_v2.py

import os
import logging
import numpy as np
import pandas as pd
from django.core.management.base import BaseCommand
from django.db import transaction
from sklearn.metrics.pairwise import cosine_similarity
from collections import defaultdict
from hormone_app.models import (
    GNPrediction, 
    HormoneReceptorFull,
    HormoneReceptorInfo
)

logger = logging.getLogger(__name__)

class Command(BaseCommand):
    help = '基于网络拓扑预测新的激素-受体相互作用（多对多版本）'

    def add_arguments(self, parser):
        parser.add_argument('--hormone-type', type=str, choices=['peptide', 'non_peptide', 'both'], 
                          default='both', help='激素类型')
        parser.add_argument('--top-k', type=int, default=50, help='每个激素预测的前K个受体')
        parser.add_argument('--threshold', type=float, default=0.1, help='预测分数阈值')
        parser.add_argument('--force', action='store_true', help='强制重新预测')

    def handle(self, *args, **options):
        self.stdout.write("🔍 开始基于网络拓扑预测激素-受体相互作用（多对多版本）...")
        
        hormone_type = options['hormone_type']
        top_k = options['top_k']
        threshold = options['threshold']
        force = options['force']
        
        if force:
            GNPrediction.objects.all().delete()
            self.stdout.write("已清空旧预测数据")
        
        total_new = 0
        
        if hormone_type in ['peptide', 'both']:
            new = self._predict_peptide(top_k, threshold)
            total_new += new
        
        if hormone_type in ['non_peptide', 'both']:
            new = self._predict_non_peptide(top_k, threshold)
            total_new += new
        
        # 显示统计
        total = GNPrediction.objects.count()
        high = GNPrediction.objects.filter(confidence_level='high').count()
        medium = GNPrediction.objects.filter(confidence_level='medium').count()
        low = GNPrediction.objects.filter(confidence_level='low').count()
        
        self.stdout.write(self.style.SUCCESS(
            f"\n✅ 预测完成！\n"
            f"本次新增: {total_new}\n"
            f"总预测数: {total}\n"
            f"高置信度: {high}\n"
            f"中置信度: {medium}\n"
            f"低置信度: {low}"
        ))
    
    def _predict_peptide(self, top_k, threshold):
        """肽类激素预测 - 多对多关系"""
        self.stdout.write("\n📊 处理肽类激素（多对多关系）...")
        
        # 获取所有肽类激素记录
        items = list(HormoneReceptorFull.objects.all())
        
        if not items:
            self.stdout.write(self.style.WARNING("没有找到肽类激素数据"))
            return 0
        
        self.stdout.write(f"总记录数: {len(items)}")
        
        # === 关键修改：构建多对多关系 ===
        # 一个激素可以通过多个受体起作用
        # 例如：同一个激素名称可能出现在多条记录中，对应不同受体
        
        hormone_to_receptors = defaultdict(set)
        receptor_to_hormones = defaultdict(set)
        
        # 按激素名称分组（而不是按ID）
        hormone_by_name = defaultdict(list)
        for item in items:
            if item.hormone_name:
                hormone_by_name[item.hormone_name.strip()].append(item)
        
        self.stdout.write(f"唯一激素名称数: {len(hormone_by_name)}")
        
        # 为每个激素名称收集它对应的所有受体
        hormone_names = list(hormone_by_name.keys())
        hormone_to_idx = {name: i for i, name in enumerate(hormone_names)}
        
        # 收集所有唯一受体
        all_receptors = []
        receptor_to_idx = {}
        
        for items_list in hormone_by_name.values():
            for item in items_list:
                receptor_key = f"{item.receptor_name}_{item.receptor_species_name}"
                if receptor_key not in receptor_to_idx:
                    receptor_to_idx[receptor_key] = len(all_receptors)
                    all_receptors.append({
                        'name': item.receptor_name,
                        'species': item.receptor_species_name,
                        'uniprot': item.receptor_uniprot_id,
                        'gene': item.receptor_coding_genes,
                        'id': item.id  # 使用其中一个记录的ID作为代表
                    })
        
        self.stdout.write(f"唯一受体数: {len(all_receptors)}")
        
        # 构建相互作用矩阵
        n_h = len(hormone_names)
        n_r = len(all_receptors)
        interaction_matrix = np.zeros((n_h, n_r))
        
        # 填充矩阵：如果激素名称对应的任何记录有这个受体，则为1
        for h_name, items_list in hormone_by_name.items():
            i = hormone_to_idx[h_name]
            for item in items_list:
                receptor_key = f"{item.receptor_name}_{item.receptor_species_name}"
                if receptor_key in receptor_to_idx:
                    j = receptor_to_idx[receptor_key]
                    interaction_matrix[i, j] = 1
        
        known_interactions = np.sum(interaction_matrix)
        self.stdout.write(f"已知相互作用数: {int(known_interactions)}")
        
        if known_interactions < 10:
            self.stdout.write(self.style.WARNING("已知相互作用太少"))
            return 0
        
        # 计算激素相似性（基于共享受体模式）
        hormone_similarity = cosine_similarity(interaction_matrix)
        
        # 预测新相互作用
        predictions = np.dot(hormone_similarity, interaction_matrix)
        
        # 过滤已知相互作用
        for i in range(n_h):
            for j in range(n_r):
                if interaction_matrix[i, j] == 1:
                    predictions[i, j] = 0
        
        # 生成预测
        created_count = 0
        with transaction.atomic():
            for i, h_name in enumerate(hormone_names):
                h_predictions = predictions[i]
                top_indices = np.argsort(h_predictions)[-top_k:][::-1]
                
                for j in top_indices:
                    score = h_predictions[j]
                    if score < threshold:
                        continue
                    
                    receptor = all_receptors[j]
                    
                    # 确定置信度
                    if score > 0.5:
                        confidence = 'high'
                    elif score > 0.3:
                        confidence = 'medium'
                    else:
                        confidence = 'low'
                    
                    # 使用该激素名称下的第一个记录的ID作为代表
                    first_item = hormone_by_name[h_name][0]
                    
                    # 检查是否已存在
                    exists = GNPrediction.objects.filter(
                        hormone_peptide_id=first_item.id,
                        hormone_type='peptide',
                        receptor_peptide_id=receptor['id'],
                        receptor_type='peptide'
                    ).exists()
                    
                    if not exists:
                        GNPrediction.objects.create(
                            hormone_peptide_id=first_item.id,
                            hormone_type='peptide',
                            receptor_peptide_id=receptor['id'],
                            receptor_type='peptide',
                            prediction_score=float(score),
                            confidence_level=confidence,
                            model_version='network_multitov1.0'
                        )
                        created_count += 1
                        
                        if created_count % 100 == 0:
                            self.stdout.write(f"  已创建 {created_count} 条预测...")
        
        self.stdout.write(f"✅ 肽类新增预测: {created_count} 条")
        return created_count
    
    def _predict_non_peptide(self, top_k, threshold):
        """非肽类激素预测 - 多对多关系"""
        self.stdout.write("\n📊 处理非肽类激素（多对多关系）...")
        
        items = list(HormoneReceptorInfo.objects.all())
        
        if not items:
            self.stdout.write(self.style.WARNING("没有找到非肽类激素数据"))
            return 0
        
        self.stdout.write(f"总记录数: {len(items)}")
        
        # 按激素名称分组
        hormone_by_name = defaultdict(list)
        for item in items:
            if item.hormone_name:
                hormone_by_name[item.hormone_name.strip()].append(item)
        
        self.stdout.write(f"唯一激素名称数: {len(hormone_by_name)}")
        
        hormone_names = list(hormone_by_name.keys())
        hormone_to_idx = {name: i for i, name in enumerate(hormone_names)}
        
        # 收集所有唯一受体
        all_receptors = []
        receptor_to_idx = {}
        
        for items_list in hormone_by_name.values():
            for item in items_list:
                receptor_key = f"{item.receptor_name}_{item.receptor_species_name}"
                if receptor_key not in receptor_to_idx:
                    receptor_to_idx[receptor_key] = len(all_receptors)
                    all_receptors.append({
                        'name': item.receptor_name,
                        'species': item.receptor_species_name,
                        'uniprot': item.receptor_uniprot_id,
                        'gene': item.receptor_coding_genes,
                        'id': item.id
                    })
        
        self.stdout.write(f"唯一受体数: {len(all_receptors)}")
        
        # 构建相互作用矩阵
        n_h = len(hormone_names)
        n_r = len(all_receptors)
        interaction_matrix = np.zeros((n_h, n_r))
        
        for h_name, items_list in hormone_by_name.items():
            i = hormone_to_idx[h_name]
            for item in items_list:
                receptor_key = f"{item.receptor_name}_{item.receptor_species_name}"
                if receptor_key in receptor_to_idx:
                    j = receptor_to_idx[receptor_key]
                    interaction_matrix[i, j] = 1
        
        known_interactions = np.sum(interaction_matrix)
        self.stdout.write(f"已知相互作用数: {int(known_interactions)}")
        
        if known_interactions < 10:
            self.stdout.write(self.style.WARNING("已知相互作用太少"))
            return 0
        
        # 计算激素相似性
        hormone_similarity = cosine_similarity(interaction_matrix)
        
        # 预测
        predictions = np.dot(hormone_similarity, interaction_matrix)
        
        # 过滤已知
        for i in range(n_h):
            for j in range(n_r):
                if interaction_matrix[i, j] == 1:
                    predictions[i, j] = 0
        
        # 生成预测
        created_count = 0
        with transaction.atomic():
            for i, h_name in enumerate(hormone_names):
                h_predictions = predictions[i]
                top_indices = np.argsort(h_predictions)[-top_k:][::-1]
                
                for j in top_indices:
                    score = h_predictions[j]
                    if score < threshold:
                        continue
                    
                    receptor = all_receptors[j]
                    
                    if score > 0.5:
                        confidence = 'high'
                    elif score > 0.3:
                        confidence = 'medium'
                    else:
                        confidence = 'low'
                    
                    first_item = hormone_by_name[h_name][0]
                    
                    exists = GNPrediction.objects.filter(
                        hormone_non_peptide_id=first_item.id,
                        hormone_type='non_peptide',
                        receptor_non_peptide_id=receptor['id'],
                        receptor_type='non_peptide'
                    ).exists()
                    
                    if not exists:
                        GNPrediction.objects.create(
                            hormone_non_peptide_id=first_item.id,
                            hormone_type='non_peptide',
                            receptor_non_peptide_id=receptor['id'],
                            receptor_type='non_peptide',
                            prediction_score=float(score),
                            confidence_level=confidence,
                            model_version='network_multitov1.0'
                        )
                        created_count += 1
        
        self.stdout.write(f"✅ 非肽类新增预测: {created_count} 条")
        return created_count