# 和输出unknown的原版代码相比，我们此版代码只修改corboundary_print函数，进而不输出包含unknown的结果
import os
import sys
import re
import logging
import multiprocessing
from multiprocessing import Pool
from collections import defaultdict, Counter

import pysam
from . import env
from .align import *
from .logger import ProgressBar
from .utils import grouper,tree

LOGGER = logging.getLogger('ECCTRY')

import itertools
import operator
import collections
import pyfaidx
import copy

def are_strands_equal(tuple1, tuple2):
    # 如果长度不同，直接返回 False
    if len(tuple1) != len(tuple2):
        return False

    # 将元组转换为列表，方便操作
    list1 = list(tuple1)
    list2 = list(tuple2)
    
    # 生成 list1 的负值版本
    neg_list1 = [-x for x in list1]
    
    # 检查四种情况：
    # 1. list1 与 list2
    # 2. list1 与 reversed_list2
    # 3. neg_list1 与 list2
    # 4. neg_list1 与 reversed_list2
    
    reversed_list2 = list2[::-1]
    
    for lst in [list1, neg_list1]:
        # 复制列表用于轮换操作
        current = lst.copy()
        
        # 检查所有可能的轮换
        for _ in range(len(lst)):
            if current == list2 or current == reversed_list2:
                return True
            # 轮换一次
            current.insert(0, current.pop())
    
    return False


def are_intervals_equal(struct1, struct2,thre = 0.8):
    # 判断intervals_tuple是否能聚到一起
    if len(struct1) != len(struct2):
        return False

    list_struct1 = list(struct1)
    list_struct2 = list(struct2)

    intervals1 = [(list_struct1[i], list_struct1[i+1]) for i in range(0, len(list_struct1), 2)]
    intervals2 = [(list_struct2[i], list_struct2[i+1]) for i in range(0, len(list_struct2), 2)]

    for (s1, e1), (s2, e2) in zip(intervals1, intervals2):

        len1 = abs(e1 - s1)
        len2 = abs(e2 - s2)

        ###[Q]为什么会出现长度为0的情况，还没有进行相关检查
        if len1 == 0 or len2 == 0:
            return False

        overlap_start = max(s1, s2)
        overlap_end = min(e1, e2)
        overlap_len = max(0, overlap_end - overlap_start)

        min_interval_len = min(len1, len2)  
        overlap_ratio = overlap_len / min_interval_len  

        if overlap_ratio < thre:
            return False
    
    # 所有区间都满足条件
    return True

def parse_chr_name(chr_name):
    """解析染色体名称，处理数字、X、Y和其他特殊名称"""
    # 移除"chr"前缀
    chr_name = chr_name.replace("chr", "")
    
    # 处理数字染色体
    if chr_name.isdigit():
        return int(chr_name)
    
    # 处理X、Y染色体
    if chr_name in ["X", "Y"]:
        return 23 if chr_name == "X" else 24
    
    # 处理其他特殊染色体（如KI270733.1）
    # 为每个唯一的特殊染色体分配一个唯一的整数值
    # 使用一个全局字典来跟踪已分配的特殊染色体
    global special_chr_map
    if not hasattr(parse_chr_name, 'special_chr_map'):
        parse_chr_name.special_chr_map = {}
    
    if chr_name not in parse_chr_name.special_chr_map:
        # 从25开始为特殊染色体分配ID
        next_id = len(parse_chr_name.special_chr_map) + 25
        parse_chr_name.special_chr_map[chr_name] = next_id
    
    return parse_chr_name.special_chr_map[chr_name]

def keywork_sort(combined_keys,intervals_list):

    if len(intervals_list) != 2 * len(combined_keys):
        raise ValueError(f"intervals_list长度必须是combined_keys的2倍")

    elements = []
    for idx in range(len(combined_keys)):
        combined_val = combined_keys[idx]
        # 提取当前combined对应的intervals start位置（intervals_list中第2*idx位）
        interval_start = intervals_list[2 * idx]
        # 排序关键字：(绝对值，start位置)，用于确定基准元素
        sort_key = (abs(combined_val), interval_start)
        # 存储：(排序关键字, combined值, 对应intervals的[start, end])
        elements.append((sort_key, combined_val, intervals_list[2*idx : 2*(idx+1)]))

    elements_sorted_by_key = sorted(elements, key=lambda x: x[0])
    base_element = elements_sorted_by_key[0]
    base_idx = combined_keys.index(base_element[1])
    sorted_combined = combined_keys[base_idx:] + combined_keys[:base_idx]
    sorted_intervals = intervals_list[2*base_idx:] + intervals_list[:2*base_idx]

    if sorted_combined[0] >= 0:
        return sorted_combined, sorted_intervals
    else:
        negated_combined = [-x for x in sorted_combined]
        reversed_combined = negated_combined[::-1]
        final_combined = [reversed_combined[-1]] + reversed_combined[:-1]

        interval_pairs = [sorted_intervals[i:i+2] for i in range(0,len(sorted_intervals),2)]
        reversed_interval_pairs = interval_pairs[::-1]
        final_interval_pairs = [reversed_interval_pairs[-1]] + reversed_interval_pairs[:-1]
        final_intervals = [val for pair in final_interval_pairs for val in pair]

    return final_combined, final_intervals


def get_interval_list(final_combined, final_intervals):

    current_interval_list = []

    for idx in range(len(final_combined)):
        combined_val = final_combined[idx]
        start = final_intervals[2*idx]      # final_intervals取第2*idx个元素（start）
        end = final_intervals[2*idx + 1]    # final_intervals取第2*idx+1个元素（end）
        sorted_group = [combined_val, start, end]  # 有序组：1+2元素组合
        current_interval_list.append(sorted_group)

    return current_interval_list


def cluster_chunk(all_ret_chunk):
    clustered_tree = tree()  # 假设tree()返回嵌套字典结构（如defaultdict）
    
    for read_info in all_ret_chunk:
        read_id, segments, ccs, raw, fraggroup = read_info

        # 20251215，出现IndexError: list index out of range，base_element = elements_sorted_by_key[0]进行的修改
        if not fraggroup.intervals:
            continue

        # 生成关键值
        combined_keys = []
        intervals_list = []

        for interval in fraggroup.intervals:
            ctg_num = parse_chr_name(interval.chr)
            combined = ctg_num * interval.strand
            combined_keys.append(combined)

            intervals_list.append(interval.start)
            intervals_list.append(interval.end)

        final_combined, final_intervals = keywork_sort(combined_keys,intervals_list)

        combined_tuple = tuple(final_combined)
        intervals_tuple = tuple(final_intervals)

        current_interval_list = get_interval_list(final_combined, final_intervals)


        if combined_tuple in clustered_tree:
            matched = False
            for existing_interval in clustered_tree[combined_tuple]:
                if are_intervals_equal(intervals_tuple, existing_interval):
                    # 追加键值并且追加名字 
                    existing_interval_list, existing_read_ids = clustered_tree[combined_tuple][existing_interval]
                    for idx, curr_interval in enumerate(current_interval_list):
                        existing_interval_list[idx].append(curr_interval)
                    existing_read_ids.append(read_id)
                    matched = True
                    break
            if not matched:
                init_interval_list = [[curr] for curr in current_interval_list]
                init_read_ids = [read_id]
                clustered_tree[combined_tuple][intervals_tuple] = [init_interval_list, init_read_ids]
        else:
            init_interval_list = [[curr] for curr in current_interval_list]
            init_read_ids = [read_id]
            clustered_tree[combined_tuple] = {intervals_tuple: [init_interval_list, init_read_ids]}


    return clustered_tree

def merge_interval_list(c2_dict, c1, i1):
    match_indices = []

    for c1_idx, c1_chr in enumerate(c1):
        i1_start_idx = 2 * c1_idx
        i1_end_idx = i1_start_idx + 2
        i1_interval = i1[i1_start_idx:i1_end_idx]  # 格式：(star, end) 

        if c1_chr not in c2_dict:
            return False, None            

        found = False
        for c2_item in c2_dict[c1_chr]:
            c2_interval = c2_item[:2]
            c2_idx = c2_item[2] 

            if c2_idx in match_indices:
                continue  # 跳过重复索引，继续寻找下一个匹配

            if are_intervals_equal(i1_interval, c2_interval):
                match_indices.append(c2_idx)
                found = True
                break

        if not found:
            return False, None

    return True, match_indices


# 合并子集分支到父集分支
def merge_subset_tree(clustered_tree):

    from collections import defaultdict

    tree_items = []   
    for combined_tuple, interval_dict in clustered_tree.items():
        for intervals_tuple, (interval_list, read_ids) in interval_dict.items():
            tree_items.append({
                "combined": combined_tuple,
                "intervals": intervals_tuple,
                "interval_list": [sublist.copy() for sublist in interval_list],  # 深拷贝避免修改原数据
                "read_ids": read_ids.copy()
            })

    merged_tree = tree()
    for item in tree_items:
        combined = item["combined"]
        intervals = item["intervals"]
        merged_tree[combined][intervals] = [item["interval_list"], item["read_ids"]]

    to_remove = set()

    # 双重遍历所有分支对（i=子集候选，j=父集候选）
    for i in range(len(tree_items)):
        item1 = tree_items[i]
        c1 = item1["combined"]  # 子集 combined_tuple
        i1 = item1["intervals"]  # 子集 intervals_tuple
        il1 = item1["interval_list"]  # 子集 interval_list
        rid1 = item1["read_ids"]  # 子集 read_id_list

        if (c1, i1) in to_remove:
            continue  

        for j in range(len(tree_items)):
            if i == j:
                continue  # 跳过自身或已标记删除的父集

            item2 = tree_items[j]
            c2 = item2["combined"]  # 父集 combined_tuple
            i2 = item2["intervals"]  # 父集 intervals_tuple
            il2 = item2["interval_list"]  # 父集 interval_list
            rid2 = item2["read_ids"]  # 父集 read_id_list

            # 条件1：父集长度必须≥子集长度（否则无法包含）
            if (c2, i2) in to_remove or len(c2) < len(c1):
                continue

            # 条件2：按顺序匹配 c1 的所有元素在 c2 中的索引（支持重复元素）
            c2_dict = {}
            for idx, chr_val in enumerate(c2):
                start_idx = 2 * idx
                end_idx = 2 *idx + 2
                interval = i2[start_idx:end_idx] + (idx,)
                if chr_val in c2_dict:
                    c2_dict[chr_val].append(interval)
                else:
                    c2_dict[chr_val] = [interval]

            ###############ing###############################

            merge_situ, merge_idx = merge_interval_list(c2_dict, c1, i1)
            if not merge_situ:
                c1_neg = tuple(-chr_val for chr_val in c1)
                merge_situ, merge_idx = merge_interval_list(c2_dict, c1_neg, i1)

            if merge_situ:
                merged_il2 = [sublist.copy() for sublist in il2]

                for indx, idx in enumerate(merge_idx):
                    merged_il2[idx].extend(il1[indx])

                merged_rid = rid2.copy()
                merged_rid.extend(rid1)

                merged_tree[c2][i2] = [merged_il2, merged_rid]

                to_remove.add((c1, i1))

                break

    # 删除所有标记为待删除的子分支
    for (combined, intervals) in to_remove:
        if combined in merged_tree and intervals in merged_tree[combined]:
            del merged_tree[combined][intervals]
            # 若该combined下已无intervals，删除整个combined键
            if not merged_tree[combined]:
                del merged_tree[combined]

    return merged_tree


def get_mode(values, default=0, is_start=True):
    """
    计算众数，支持多个众数时的特殊处理：
    :param values: 待计算众数的列表
    :param default: 无数据时的默认值
    :param is_start: 是否为开始位置（True→多个众数取最小；False→多个众数取最大）
    :return: 最终众数
    """
    if not values:
        return default
    # 统计每个值的频次
    value_counts = Counter(values)
    # 找到最大频次
    max_freq = max(value_counts.values())
    # 收集所有频次等于最大频次的众数
    modes = [val for val, freq in value_counts.items() if freq == max_freq]
    # 多个众数时的处理逻辑
    if len(modes) > 1:
        if is_start:
            return min(modes)  # 开始位置：取最小众数
        else:
            return max(modes)  # 结束位置：取最大众数
    # 单个众数时直接返回
    return modes[0]


def corboundary_print(merged_tree, out_dir, prefix):

    output_file = f"{out_dir}/{prefix}_eccdna_information.txt"
    

    chr_name_map = {
        1: "chr1", 2: "chr2", 3: "chr3", 4: "chr4", 5: "chr5", 
        6: "chr6", 7: "chr7", 8: "chr8", 9: "chr9", 10: "chr10",
        11: "chr11", 12: "chr12", 13: "chr13", 14: "chr14", 15: "chr15",
        16: "chr16", 17: "chr17", 18: "chr18", 19: "chr19", 20: "chr20",
        21: "chr21", 22: "chr22", 23: "chrX", 24: "chrY"
    }

    counter = 1
    filtered_count = 0

    with open(output_file, 'w', encoding='utf-8') as f:
        # 写入表头（与参考函数一致）
        header = (
            "id\tmerge_region\tmerge_len\tnum_region\t"
            "eccdna_status\tnumreads\n"
        )
        f.write(header)


        for combined_key, intervals_dict in merged_tree.items():
            # 遍历第二层（intervals_tuple：如 (10,20,50,60,30,40)）
            for intervals_tuple, (interval_list, read_ids) in intervals_dict.items():
                # 存储当前分支的所有区域字符串（如 "chr5:123-456_+"）
                hit_strings = []
                merge_len = 0  # 所有区域长度之和
                region_num = 0  # 区域数量（即 interval_list 的元素个数）
                has_unknown_chr = False


                # 遍历 interval_list 中的每个元素A（对应一个区间）
                for sampos_list in interval_list:

                    chr_num_raw = sampos_list[0][0]  # sampos_list[0][0] 是染色体编号（整数，可能为负）
                    chr_abs_num = abs(chr_num_raw)  # 取绝对值（用于映射染色体名称）
                    strand_sym = '+' if chr_num_raw > 0 else '-'  # 正→+，负→-

                    chr_name = chr_name_map.get(chr_abs_num, None)

                    if chr_name is None:
                        has_unknown_chr = True
                        filtered_count += 1
                        LOGGER.warning(f"检测到未知染色体ID {chr_abs_num}，跳过该分支：combined_key={combined_key}, intervals_tuple={intervals_tuple}")
                        break  # 跳出当前interval_list循环，不处理该分支


                    start_values = []
                    end_values = []

                    for sublist in sampos_list:
                        start_values.append(sublist[1])  # 第2位 → start
                        end_values.append(sublist[2])    # 第3位 → end

                    start_mode = get_mode(start_values, default=0, is_start=True)
                    end_mode = get_mode(end_values, default=0, is_start=False)

                    # 步骤4：确保 start ≤ end（避免区间反转）
                    if start_mode > end_mode:
                        start_mode, end_mode = end_mode, start_mode
                        

                    # 步骤5：构建区域字符串（如 "chr5:124089936-124090338_+"）
                    hit_str = f"{chr_name}:{start_mode}-{end_mode}_{strand_sym}"
                    hit_strings.append(hit_str)

                    # 步骤6：计算当前区域长度并累加
                    region_len = end_mode - start_mode + 1  # 闭区间长度
                    merge_len += region_len
                    region_num += 1

                if has_unknown_chr:
                    continue

                # 步骤7：处理当前分支的所有信息
                # 拼接 merge_region（多个区域用逗号分隔）
                merge_region = ",".join(hit_strings) if hit_strings else "NA"
                numreads = len(read_ids)  # 读段数量
                eccdna_status = "ccs_cyclic"

                # 步骤8：构建输出行
                line = (
                    f"ec{counter}\t"
                    f"{merge_region}\t"
                    f"{merge_len}\t"
                    f"{region_num}\t"
                    f"{eccdna_status}\t"
                    f"{numreads}\n"
                )

                # 写入文件
                f.write(line)
                counter += 1  # ID自增
                
    LOGGER.info(f"输出文件已保存至：{output_file}")
    LOGGER.info(f"共处理 {counter-1} 个有效分支，过滤掉 {filtered_count} 个包含未知染色体的分支")

    return


def merge_trees(tree1, tree2, interval_thre=0.8):
    """
    合并两个eccDNA聚类树（无深拷贝，直接修改tree1）
    :param tree1: 第一个聚类树（嵌套字典结构，会被直接修改）
    :param tree2: 第二个聚类树（嵌套字典结构，原数据保留）
    :param interval_thre: 区间匹配的重叠阈值（默认0.8）
    :return: 合并后的树（即修改后的tree1）
    """
    # 直接使用tree1作为合并后的树，无深拷贝
    merged_tree = tree1

    # 遍历tree2的所有一级键（combined_tuple）
    for combined_tuple in tree2:
        # 情况1：tree1中已存在该combined_tuple → 尝试匹配intervals_tuple
        if combined_tuple in merged_tree:
            tree2_intervals_dict = tree2[combined_tuple]
            merged_intervals_dict = merged_tree[combined_tuple]

            # 遍历tree2当前combined_tuple下的所有intervals_tuple
            for intervals_tuple2, (interval_list2, read_ids2) in tree2_intervals_dict.items():
                matched = False
                # 遍历merged_tree中当前combined_tuple下的intervals_tuple，尝试匹配
                for intervals_tuple1, (interval_list1, read_ids1) in merged_intervals_dict.items():
                    if are_intervals_equal(intervals_tuple1, intervals_tuple2, thre=interval_thre):
                        # 合并interval_list：直接追加tree2的列表（零拷贝）
                        if len(interval_list1) != len(interval_list2):
                            continue
                        for idx in range(len(interval_list1)):
                            # 直接extend tree2的原列表，无任何拷贝
                            interval_list1[idx].extend(interval_list2[idx])
                        
                        # 合并read_ids：直接去重追加tree2的原列表（零拷贝）
                        for rid in read_ids2:
                                read_ids1.append(rid)
                        
                        matched = True
                        break
                
                # 情况1.2：无匹配的intervals_tuple → 直接复用tree2的二级键数据（零拷贝）
                if not matched:
                    # 直接将tree2的列表引用赋值给merged_tree，无拷贝
                    merged_intervals_dict[intervals_tuple2] = [interval_list2, read_ids2]
        
        # 情况2：tree1中不存在该combined_tuple → 直接复用tree2的一级键数据（零拷贝）
        else:
            # 直接将tree2的intervals_dict引用赋值给merged_tree，无拷贝
            merged_tree[combined_tuple] = tree2[combined_tuple]

    return merged_tree


def cluster_all_ret(all_ret, threads):
    # 分割数据
    chunks = []
    total_items = len(all_ret)
    
    # 处理线程数大于数据量的情况
    if threads > total_items:
        threads = total_items
        chunk_size = 1
        remainder = 0
    else:
        chunk_size = total_items // threads
        remainder = total_items % threads
    
    start_idx = 0
    for i in range(threads):
        # 前remainder个块多分配1个元素
        end_idx = start_idx + chunk_size + (1 if i < remainder else 0)
        chunks.append(all_ret[start_idx:end_idx])
        start_idx = end_idx

    # 创建进程池
    with Pool(processes=threads) as pool:
        # 并行处理每个块（不需要partial，因为cluster_chunk只需要一个参数）
        results = pool.map(cluster_chunk, chunks)

    # 合并所有局部的聚类树
    final_tree = tree()
    for result in results:
        final_tree = merge_trees(final_tree, result)

    return final_tree


def print_tree(tree, out_dir, prefix, max_interval_list_len=500):
    """
    格式化打印eccDNA聚类树到指定文件，按列输出核心信息
    :param tree: 待打印的聚类树（嵌套字典结构）
    :param out_dir: 输出目录（自动创建不存在的目录）
    :param prefix: 输出文件前缀（最终文件名：{prefix}_tree.txt）
    :param max_interval_list_len: interval_list字符串最大长度（避免输出过长，默认500）
    """
    # 确保输出目录存在
    os.makedirs(out_dir, exist_ok=True)
    # 定义输出文件路径
    output_file = f"{out_dir}/{prefix}_tree.txt"

    # 打开文件并写入内容（utf-8编码，覆盖写入）
    with open(output_file, 'w', encoding='utf-8') as f:
        # 1. 写入表头（制表符分隔）
        header = "\t".join([
            "combined_tuple",
            "intervals_tuple",
            "interval_list",
            "read_ids"
        ])
        f.write(header + "\n")
        f.write("-" * 120 + "\n")  # 写入分隔线

        # 2. 遍历树结构，逐行写入数据
        for combined_tuple in tree:
            intervals_dict = tree[combined_tuple]
            # 处理空intervals_dict的情况
            if not intervals_dict:
                row = "\t".join([
                    str(combined_tuple),
                    "无",
                    "无",
                    "无"
                ])
                f.write(row + "\n")
                continue

            # 遍历每个intervals_tuple
            for intervals_tuple, (interval_list, read_ids) in intervals_dict.items():
                # 处理interval_list：格式化+超长截断
                interval_list_str = str(interval_list)
                if len(interval_list_str) > max_interval_list_len:
                    interval_list_str = interval_list_str[:max_interval_list_len] + "..."
                
                # 处理read_ids：列表转逗号分隔字符串
                read_ids_str = ",".join(read_ids) if isinstance(read_ids, list) else str(read_ids)
                
                # 拼接行数据（制表符分隔）并写入
                row = "\t".join([
                    str(combined_tuple),
                    str(intervals_tuple),
                    interval_list_str,
                    read_ids_str
                ])
                f.write(row + "\n")

        # 3. 写入统计信息
        total_combined = len(tree)
        total_intervals = sum(len(intervals_dict) for intervals_dict in tree.values())
        f.write("-" * 120 + "\n")
        f.write(f"统计信息：共 {total_combined} 个combined_tuple分支，{total_intervals} 个intervals_tuple子分支\n")

# def print_cluster_tree(final_tree, out_dir, prefix):
#     """
#     打印聚类树内容，每行格式：ctg_strand_tuple \t intervals_tuple \t read_id \t ctg_tuple \t positions_tuple \n
#     适配fraggroup结构，从FragmentGroup中提取染色体和位置信息
#     """
#     # 输出文件路径（与其他结果文件放在同一目录）
#     output_file = f"{out_dir}/tmp/{prefix}_tree.txt"

#     with open(output_file, 'w') as f:
#         # 遍历两层嵌套结构（combined_tuple → intervals_tuple → read_info_list）
#         for combined_tuple in final_tree:
#             for intervals_tuple in final_tree[combined_tuple]:
#                 # 遍历该分组下的所有read_info
#                 for read_info in final_tree[combined_tuple][intervals_tuple]:
#                     read_id = read_info[0]  # read_info格式为(read_id, segments, ccs, raw, fraggroup)
#                     fraggroup = read_info[4]  # 获取fraggroup（替代原cand_hit_list）

#                     # 从fraggroup的intervals中提取染色体编号和位置信息
#                     ctg_list = []
#                     positions_list = []
#                     # 确保fraggroup和intervals有效
#                     for interval in fraggroup.intervals:

#                         strand = interval.strand

#                         # 解析染色体编号（使用parse_chr_name转换为数字）
#                         ctg_num = parse_chr_name(interval.chr)
#                         ctg_list.append(ctg_num * strand)
                        
#                         # 提取起始和结束位置（使用FragmentGroup的interval属性）
#                         positions_list.append(interval.start)
#                         positions_list.append(interval.end)

#                     # 转换为元组并格式化
#                     ctg_tuple = str(tuple(ctg_list))
#                     positions_tuple = str(tuple(positions_list))

#                     # 格式化元组（便于阅读）
#                     formatted_combined = str(combined_tuple)
#                     formatted_intervals = str(intervals_tuple)
                    
#                     # 按格式写入一行
#                     line = f"{formatted_combined}\t{formatted_intervals}\t{read_id}\t{ctg_tuple}\t{positions_tuple}\n"
#                     f.write(line)

# def print_cluster_tree_outputformat(clustered_tree, out_dir, prefix):
#     """
#     打印聚类树内容，每行格式：
#     id	merge_region	merge_len	num_region	eccdna_status	numreads	totalbase	coverage	consensus_len	consensus_status	read_name
#     ec1	chr5:124089936-124090338_+,chr4:148751511-148753435_+	2328	2	ccs_cyclic	5	0	0	0	0	GN234_1,GN235_2...
#     """
#     output_file = f"{out_dir}/tmp/{prefix}_ccs_candidate_ecc.txt"

#     # 定义染色体名称正向映射（数字→染色体名称）
#     chr_name_map = {
#         1: "chr1", 2: "chr2", 3: "chr3", 4: "chr4", 5: "chr5", 
#         6: "chr6", 7: "chr7", 8: "chr8", 9: "chr9", 10: "chr10",
#         11: "chr11", 12: "chr12", 13: "chr13", 14: "chr14", 15: "chr15",
#         16: "chr16", 17: "chr17", 18: "chr18", 19: "chr19", 20: "chr20",
#         21: "chr21", 22: "chr22", 23: "chrX", 24: "chrY"
#     }
    
#     # 获取特殊染色体的反向映射（编号→原始名称）
#     special_chr_reverse_map = {}
#     if hasattr(parse_chr_name, 'special_chr_map'):
#         special_chr_reverse_map = {v: k for k, v in parse_chr_name.special_chr_map.items()}
    
#     # 初始化ID计数器
#     counter = 1
    
#     with open(output_file, 'w', encoding='utf-8') as f:
#         header = (
#             "id\tmerge_region\tmerge_len\tnum_region\t"
#             "eccdna_status\tnumreads\ttotalbase\tcoverage\t"
#             "consensus_len\tconsensus_status\tread_name\n"
#         )
#         f.write(header)
#         # 遍历聚类树的一级节点（combined_key）
#         for combined_key, intervals_dict in clustered_tree.items():
#             # 遍历二级节点（intervals_tuple）
#             for intervals_tuple, read_info_list in intervals_dict.items():
#                 # 使用第一个read的fraggroup作为代表
#                 read_info_first = read_info_list[0]
#                 read_id, segments, ccs, raw, fraggroup = read_info_first  # 替换为fraggroup

#                 hit_strings = []
#                 merge_len = 0
#                 region_num = 0

#                 # 从fraggroup的intervals提取信息（替代原cand_hit_list）
#                 for interval in fraggroup.intervals:
#                     # 确定链方向符号
#                     strand_sym = '+' if interval.strand==1 else '-'
                    
#                     # 构建区域字符串（染色体:起始-结束_链方向）
#                     hit_str = f"{interval.chr}:{interval.start}-{interval.end}_{strand_sym}"
#                     hit_strings.append(hit_str)
                    
#                     # 计算区域长度并累加
#                     region_len = interval.end - interval.start + 1
#                     merge_len += region_len
                    
#                     region_num += 1

#                 # 拼接merge_region字符串
#                 merge_region = ",".join(hit_strings)

#                 # 处理读段ID列表
#                 read_ids = [read_info[0] for read_info in read_info_list]
#                 read_ids_str = ",".join(read_ids)
#                 numreads = len(read_ids)

#                 # 保留原始默认值（可根据实际需求从fraggroup补充）
#                 totalbase = 0
#                 coverage = 0
#                 consensus_len = 0
#                 consensus_status = 0
                
#                 # 构建输出行
#                 line = (
#                     f"ec{counter}\t"
#                     f"{merge_region}\t"
#                     f"{merge_len}\t"
#                     f"{region_num}\t"
#                     f"ccs_cyclic\t"
#                     f"{numreads}\t"
#                     f"{totalbase}\t"
#                     f"{coverage}\t"
#                     f"{consensus_len}\t"
#                     f"{consensus_status}\t"
#                     f"{read_ids_str}\n"
#                 )

#                 f.write(line)
#                 counter += 1  # 自增ID

# def print_fill_gap_filtered(all_fill_gap_filtered, out_dir, prefix):
#     """
#     一些比对结果包含gap，我们将对这样的结果进行调整，下面打印的是调整失败的结果
    
#     每行格式：


#     """
#     out_file = f"{out_dir}/tmp/{prefix}_fill_gap_filtered.txt"

#     with open(out_file, 'w') as f:
#         # 写入表头
#         header = "read_id\tq_st\tq_en\tctg\tstrand\tr_st\tr_en\tq_st\tq_en\tmlen\tblen\tcigar_string\n"
#         f.write(header)

#         # 遍历每个失败的 ccs_hit_list
#         for read_id, segments, ccs, raw, cand_hit_list in all_fill_gap_filtered:
#             for hit in cand_hit_list:
#                 # 构建每一行的内容
#                 line = f"{read_id}\t{hit.q_st}\t{hit.q_en}\t{hit.ctg}\t{hit.strand}\t{hit.r_st}\t{hit.r_en}\t{hit.mlen}\t{hit.blen}\t{hit.cigar}\n"
#                 f.write(line)

#     return
