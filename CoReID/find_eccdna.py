import os
import sys
import re
import logging
import multiprocessing
from multiprocessing import Pool
from collections import defaultdict

import pysam
from . import env
from .align import *
from .logger import ProgressBar
from .utils import grouper,tree

LOGGER = logging.getLogger('ECCTRY')

import itertools
import operator
import collections

# [Q]此过滤条件是合理的吗？
# 同时也过滤掉了sp？？？
def check_raw_hits(raw_hit_list, raw, ccs):
    """检查原始比对是否通过过滤条件"""
    for raw_hit in raw_hit_list:
        # 说明是线性的或者序列的拷贝数不够
        if raw_hit and raw_hit.mlen > max(len(raw) * 0.8, len(raw) - 200):
            return False
        ### 如果某一条回帖记录满足以下条件，说明CCS有误或者比对出现了错误
        if raw_hit and raw_hit.mlen > 1.5 * len(ccs):
            return False
    return True


def check_overlap(hit_list):
    """
    判断raw_hit_list中是否存在：
    1. 同一ctg（contig）
    2. 相同strand（链方向）
    3. [r_st, r_en]区间重叠
    的subhit，满足则返回True，否则返回False
    """

    groups = {}
    for subhit in hit_list:
        # 用(ctg, strand)作为分组键
        key = (subhit.ctg, subhit.strand)
        if key not in groups:
            groups[key] = []
        groups[key].append(subhit)

    # 对每个分组（同ctg+同strand）检查区间重叠
    for key, hits in groups.items():
        # 按r_st排序当前分组内的区间
        sorted_hits = sorted(hits, key=lambda x: x.r_st)
        # 检查相邻区间是否重叠
        for i in range(1, len(sorted_hits)):
            prev = sorted_hits[i-1]
            curr = sorted_hits[i]
            # 区间重叠条件：前一个的结束 > 当前的开始
            if prev.r_en > curr.r_st:
                return True  # 发现重叠直接返回True

    return False


def process_ccs_hit_list(ccs_hit_list):
# 从最长的开始遍历，合理吗
    if not ccs_hit_list:
        return []

    max_index, max_mlen_subhit = max(enumerate(ccs_hit_list), key=lambda x:(x[1].mlen, -x[0]))
    result_list = [max_mlen_subhit]

    for i in range(max_index + 1, len(ccs_hit_list)): 
        overlap = False
        subhit = ccs_hit_list[i]
        
        for existing in result_list:
            if (existing.ctg == subhit.ctg and
                existing.strand == subhit.strand and
                existing.r_st <= subhit.r_en and 
                existing.r_en >= subhit.r_st):
                overlap = True
                break

        if overlap:
            break  # 发现重叠则停止遍历
        result_list.append(subhit)

    return result_list

def scan_ccs_chunk(chunk, mapq):
    reads_cnt = defaultdict(int) # 统计不同过滤截断的序列数量
    ret = [] # 存储最终通过所有过滤的环形RNA候选

    for read_id, segments, ccs, raw in chunk:
        # Filter 1 - Remove linear mapped reads
        raw_hit_list = get_primary_alignment(env.ALIGNER.map(raw))

        if not check_overlap(raw_hit_list):
            reads_cnt["raw_alignment_interval_noverlap"] += 1
            continue

        if not check_raw_hits(raw_hit_list, raw, ccs):
            reads_cnt["filtered_by_raw_alignment"] += 1
            continue


        # 从原始序列中直接得到eccDNA，如果有结果的话就不使用CCS的结果
        raw_fraggroup, read_raw = raw_get_ecc(read_id, raw, raw_hit_list, mapq)

        if raw_fraggroup is not None:
            ret.append((read_id, segments, ccs, raw, raw_fraggroup.to_dict()))
            continue

        else:
            ccs_hit_list = get_primary_alignment(env.ALIGNER.map(ccs * 2))

            cand_hit_list = process_ccs_hit_list(ccs_hit_list)

            if not cand_hit_list:
                continue

            # first_q_st = cand_hit_list[0].q_st
            # last_q_en = cand_hit_list[-1].q_en
            # ccs_consensus = (ccs * 2)[first_q_st : last_q_en + 1]

            ccs_fraggroup = ccs_hit_to_fraggroup(segments, cand_hit_list)
            ret.append((read_id, segments, ccs, raw, ccs_fraggroup.to_dict()))  

    return reads_cnt, ret





def scan_ccs_reads(ccs_seq, ref_fasta, out_dir, prefix, threads, mapq):
    import mappy as mp

    faidx = Faidx(ref_fasta)
    contig_len = faidx.contig_len
    faidx.close()

    # First scanning using minimap2
    minimap_aligner = mp.Aligner(ref_fasta, n_threads=threads, preset='map-ont')

    chunk_size = 250
    jobs = []

    pool = Pool(threads, env.initializer,
                (minimap_aligner, contig_len, ref_fasta))


    for reads in grouper(list(ccs_seq), chunk_size):
        chunk = [[i, ] + ccs_seq[i] for i in reads if i is not None]
        jobs.append(pool.apply_async(scan_ccs_chunk, (chunk, mapq) ))
    pool.close()

    prog = ProgressBar()
    finished_cnt = 0

    reads_count = defaultdict(int)

    all_ret = []
    for job in jobs:
        finished_cnt += 1

        tmp_cnt, ret = job.get()
        for key, value in tmp_cnt.items():
            reads_count[key] += value

        all_ret.extend(ret)
        prog.update(100 * finished_cnt / len(jobs))


    pool.join()
    # 主进程恢复 FragmentGroup
    all_ret = [
        (r[0], r[1], r[2], r[3], FragmentGroup.from_dict(r[4]))
        for r in all_ret
    ]

    prog.update(100)

    return reads_count, all_ret


# def are_intervals_equal(struct1, struct2):
#     # 判断intervals_tuple是否能聚到一起
#     if len(struct1) != len(struct2):
#         return False
#     struct1 = list(struct1)
#     struct2 = list(struct2)
#     sorted_struct1 = sorted(struct1)
#     sorted_struct2 = sorted(struct2)

#     for elem1, elem2 in zip(sorted_struct1, sorted_struct2):
#         if abs(elem1 - elem2) > 10:
#             return False

#     return True

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


def are_intervals_equal(struct1, struct2,thre = 0.85):
    # 判断intervals_tuple是否能聚到一起
    if len(struct1) != len(struct2):
        return False

    sorted_struct1 = sorted(list(struct1))
    sorted_struct2 = sorted(list(struct2))

    intervals1 = [(sorted_struct1[i], sorted_struct1[i+1]) for i in range(0, len(sorted_struct1), 2)]
    intervals2 = [(sorted_struct2[i], sorted_struct2[i+1]) for i in range(0, len(sorted_struct2), 2)]

    for (s1, e1), (s2, e2) in zip(intervals1, intervals2):


        len1 = abs(e1 - s1)
        len2 = abs(e2 - s2)

        overlap_start = max(s1, s2)
        overlap_end = min(e1, e2)
        overlap_len = max(0, overlap_end - overlap_start)

        max_len = min(len1, len2)
        overlap_ratio = overlap_len / max_len

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


def cluster_chunk(all_ret_chunk):
    clustered_tree = tree()  # 假设tree()返回嵌套字典结构（如defaultdict）
    
    for read_info in all_ret_chunk:
        read_id, segments, ccs, raw, fraggroup = read_info

        # 生成复合查询值
        combined_keys = []
        for interval in fraggroup.intervals:
            # 解析染色体名称为数字
            ctg_num = parse_chr_name(interval.chr)
            # 结合链方向生成复合键（与原逻辑保持一致）
            strand = interval.strand
            combined = ctg_num * strand
            combined_keys.append(combined)
        combined_tuple = tuple(combined_keys)

        # 从fraggroup的intervals中提取区间信息（替代原有的hit.r_st和hit.r_en）
        intervals_list = []
        for interval in fraggroup.intervals:
            intervals_list.append(interval.start)
            intervals_list.append(interval.end)
        intervals_tuple = tuple(intervals_list)

        # 查找匹配的combined_tuple（使用are_intervals_equal判断）
        matched_combined_key = None
        for existing_combined in clustered_tree:
            if are_strands_equal(combined_tuple, existing_combined):
                matched_combined_key = existing_combined
                break

        if matched_combined_key:
            # 找到匹配的combined_key，合并聚类
            # 在匹配的combined_key对应的组中查找匹配的区间结构
            matched = False
            for existing_interval in clustered_tree[matched_combined_key]:
                if are_intervals_equal(intervals_tuple, existing_interval):
                    # 结构匹配，添加到现有组
                    clustered_tree[matched_combined_key][existing_interval].append(read_info)
                    matched = True
                    break
            if not matched:
                # 结构不匹配，创建新组
                clustered_tree[matched_combined_key][intervals_tuple] = [read_info]
        else:
            # 没有匹配的combined_key，创建新的一级节点
            clustered_tree[combined_tuple] = {intervals_tuple: [read_info]}
    
    return clustered_tree



def merge_trees(tree1, tree2):
    """合并两个聚类树，处理可能的结构匹配（基于are_intervals_equal）"""
    for combined_tuple, intervals_dict in tree2.items():
        # 查找匹配的combined_tuple（使用are_intervals_equal判断）
        matched_combined_key = None
        for existing_combined in tree1:
            if are_strands_equal(combined_tuple, existing_combined):
                matched_combined_key = existing_combined
                break

        if matched_combined_key:
            # 合并相同或相似的combined_key下的区间字典
            for intervals_tuple, read_info_list in intervals_dict.items():
                # 检查当前区间是否有匹配的已有区间
                matched_interval = None
                for existing_interval in tree1[matched_combined_key]:
                    if are_intervals_equal(intervals_tuple, existing_interval):
                        matched_interval = existing_interval
                        break
                
                if matched_interval:
                    tree1[matched_combined_key][matched_interval].extend(read_info_list)
                else:
                    tree1[matched_combined_key][intervals_tuple] = read_info_list
        else:
            # 没有匹配的combined_key，直接添加
            tree1[combined_tuple] = intervals_dict
    
    return tree1


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

# def print_cluster_tree(final_tree, out_dir, prefix):
#     """
#     打印聚类树内容，每行格式：ctg_strand_tuple \t intervals_tuple \t read_id \n
#     """
#     # 输出文件路径（与其他结果文件放在同一目录）
#     output_file = f"{out_dir}/tmp/{prefix}_tree.txt"

#     with open(output_file, 'w') as f:
#         # 遍历两层嵌套结构（combined_tuple → intervals_tuple → read_info_list）
#         for combined_tuple in final_tree:
#             for intervals_tuple in final_tree[combined_tuple]:
#                 # 遍历该分组下的所有read_info
#                 for read_info in final_tree[combined_tuple][intervals_tuple]:
#                     read_id = read_info[0]  # read_info格式为(read_id, segments, ccs, raw, cand_hit_list)
                    
#                     # 格式化元组（便于阅读）
#                     formatted_combined = str(combined_tuple)
#                     formatted_intervals = str(intervals_tuple)
                    
#                     # 按格式写入一行
#                     line = f"{formatted_combined}\t{formatted_intervals}\t{read_id}\n"
#                     f.write(line)

def print_cluster_tree(final_tree, out_dir, prefix):
    """
    打印聚类树内容，每行格式：ctg_strand_tuple \t intervals_tuple \t read_id \t ctg_tuple \t positions_tuple \n
    适配fraggroup结构，从FragmentGroup中提取染色体和位置信息
    """
    # 输出文件路径（与其他结果文件放在同一目录）
    output_file = f"{out_dir}/tmp/{prefix}_tree.txt"

    with open(output_file, 'w') as f:
        # 遍历两层嵌套结构（combined_tuple → intervals_tuple → read_info_list）
        for combined_tuple in final_tree:
            for intervals_tuple in final_tree[combined_tuple]:
                # 遍历该分组下的所有read_info
                for read_info in final_tree[combined_tuple][intervals_tuple]:
                    read_id = read_info[0]  # read_info格式为(read_id, segments, ccs, raw, fraggroup)
                    fraggroup = read_info[4]  # 获取fraggroup（替代原cand_hit_list）

                    # 从fraggroup的intervals中提取染色体编号和位置信息
                    ctg_list = []
                    positions_list = []
                    # 确保fraggroup和intervals有效
                    for interval in fraggroup.intervals:

                        strand = interval.strand

                        # 解析染色体编号（使用parse_chr_name转换为数字）
                        ctg_num = parse_chr_name(interval.chr)
                        ctg_list.append(ctg_num * strand)
                        
                        # 提取起始和结束位置（使用FragmentGroup的interval属性）
                        positions_list.append(interval.start)
                        positions_list.append(interval.end)

                    # 转换为元组并格式化
                    ctg_tuple = str(tuple(ctg_list))
                    positions_tuple = str(tuple(positions_list))

                    # 格式化元组（便于阅读）
                    formatted_combined = str(combined_tuple)
                    formatted_intervals = str(intervals_tuple)
                    
                    # 按格式写入一行
                    line = f"{formatted_combined}\t{formatted_intervals}\t{read_id}\t{ctg_tuple}\t{positions_tuple}\n"
                    f.write(line)

def print_cluster_tree_outputformat(clustered_tree, out_dir, prefix):
    """
    打印聚类树内容，每行格式：
    id	merge_region	merge_len	num_region	eccdna_status	numreads	totalbase	coverage	consensus_len	consensus_status	read_name
    ec1	chr5:124089936-124090338_+,chr4:148751511-148753435_+	2328	2	ccs_cyclic	5	0	0	0	0	GN234_1,GN235_2...
    """
    output_file = f"{out_dir}/tmp/{prefix}_ccs_candidate_ecc.txt"

    # 定义染色体名称正向映射（数字→染色体名称）
    chr_name_map = {
        1: "chr1", 2: "chr2", 3: "chr3", 4: "chr4", 5: "chr5", 
        6: "chr6", 7: "chr7", 8: "chr8", 9: "chr9", 10: "chr10",
        11: "chr11", 12: "chr12", 13: "chr13", 14: "chr14", 15: "chr15",
        16: "chr16", 17: "chr17", 18: "chr18", 19: "chr19", 20: "chr20",
        21: "chr21", 22: "chr22", 23: "chrX", 24: "chrY"
    }
    
    # 获取特殊染色体的反向映射（编号→原始名称）
    special_chr_reverse_map = {}
    if hasattr(parse_chr_name, 'special_chr_map'):
        special_chr_reverse_map = {v: k for k, v in parse_chr_name.special_chr_map.items()}
    
    # 初始化ID计数器
    counter = 1
    
    with open(output_file, 'w', encoding='utf-8') as f:
        header = (
            "id\tmerge_region\tmerge_len\tnum_region\t"
            "eccdna_status\tnumreads\ttotalbase\tcoverage\t"
            "consensus_len\tconsensus_status\tread_name\n"
        )
        f.write(header)
        # 遍历聚类树的一级节点（combined_key）
        for combined_key, intervals_dict in clustered_tree.items():
            # 遍历二级节点（intervals_tuple）
            for intervals_tuple, read_info_list in intervals_dict.items():
                # 使用第一个read的fraggroup作为代表
                read_info_first = read_info_list[0]
                read_id, segments, ccs, raw, fraggroup = read_info_first  # 替换为fraggroup

                hit_strings = []
                merge_len = 0
                region_num = 0

                # 从fraggroup的intervals提取信息（替代原cand_hit_list）
                for interval in fraggroup.intervals:
                    # 确定链方向符号
                    strand_sym = '+' if interval.strand==1 else '-'
                    
                    # 构建区域字符串（染色体:起始-结束_链方向）
                    hit_str = f"{interval.chr}:{interval.start}-{interval.end}_{strand_sym}"
                    hit_strings.append(hit_str)
                    
                    # 计算区域长度并累加
                    region_len = interval.end - interval.start + 1
                    merge_len += region_len
                    
                    region_num += 1

                # 拼接merge_region字符串
                merge_region = ",".join(hit_strings)

                # 处理读段ID列表
                read_ids = [read_info[0] for read_info in read_info_list]
                read_ids_str = ",".join(read_ids)
                numreads = len(read_ids)

                # 保留原始默认值（可根据实际需求从fraggroup补充）
                totalbase = 0
                coverage = 0
                consensus_len = 0
                consensus_status = 0
                
                # 构建输出行
                line = (
                    f"ec{counter}\t"
                    f"{merge_region}\t"
                    f"{merge_len}\t"
                    f"{region_num}\t"
                    f"ccs_cyclic\t"
                    f"{numreads}\t"
                    f"{totalbase}\t"
                    f"{coverage}\t"
                    f"{consensus_len}\t"
                    f"{consensus_status}\t"
                    f"{read_ids_str}\n"
                )

                f.write(line)
                counter += 1  # 自增ID



def print_fill_gap_filtered(all_fill_gap_filtered, out_dir, prefix):
    """
    一些比对结果包含gap，我们将对这样的结果进行调整，下面打印的是调整失败的结果
    
    每行格式：


    """
    out_file = f"{out_dir}/tmp/{prefix}_fill_gap_filtered.txt"

    with open(out_file, 'w') as f:
        # 写入表头
        header = "read_id\tq_st\tq_en\tctg\tstrand\tr_st\tr_en\tq_st\tq_en\tmlen\tblen\tcigar_string\n"
        f.write(header)

        # 遍历每个失败的 ccs_hit_list
        for read_id, segments, ccs, raw, cand_hit_list in all_fill_gap_filtered:
            for hit in cand_hit_list:
                # 构建每一行的内容
                line = f"{read_id}\t{hit.q_st}\t{hit.q_en}\t{hit.ctg}\t{hit.strand}\t{hit.r_st}\t{hit.r_en}\t{hit.mlen}\t{hit.blen}\t{hit.cigar_str}\n"
                f.write(line)

    return





######################################## 从raw中是识别eccDNA #################################################################################

def multimode(data):
    counts = collections.Counter(iter(data)).most_common()
    _, mode_items = next(itertools.groupby(counts, key=operator.itemgetter(1)), (0, []))
    return list(map(operator.itemgetter(0), mode_items))

class Interval(object):
    def __init__(self, chr, start, end, strand):
        self.chr = chr
        self.start = start
        self.end = end
        self.strand = strand


class Fragment(object):
    def __init__(self, start, end, chr, chrStart, chrEnd, strand, mapQual, cigar):
        # Both start and end are 0-based and inclusive
        self.start = start
        self.end = end
        self.interval = Interval(chr, chrStart, chrEnd, strand)
        self.mapQual = mapQual
        self.cigar = cigar


#管理一组片段，用于识别滚环复制的模式，并生成每个分组的代表性参考区域。
#FragmentGroup，整数，两个片段的 start/end 位置允许的最大偏差
class FragmentGroup(object):
    def __init__(self, maxOffset=10):
        self.group = [] # 二维列表，存储分组结果（每个子列表对应一个分组，包含该分组的所有Fragment对象）。
        self.intervals = [] # 列表，存储每个分组的代表性Interval（参考基因组区域，后续通过reduce方法生成）。
        self.nfrag = 0 # 分组数目（滚环模式的 “片段类型数”，如 “片段 A→片段 B” 的循环中nfrag=2）。
        self.passes = 0 # 整数，滚环复制的完整周期数
        self.pointer = 0 # 整数，当前分组指针（标记下一个片段应加入的分组索引）
        self.addFirst = False # 布尔值，标记首尾片段是否成功加入分组。
        self.addLast = False
        self.maxOffset = maxOffset

    def _sameLocation(self, index, interval, withStart=True, withEnd=True):
        if len(self.group[index]) == 0:
            return True
        if interval.chr != self.group[index][0].interval.chr:
            return False
        if withStart:
            start = min(multimode([frag.interval.start for frag in self.group[index]]))
            if abs(start - interval.start) > self.maxOffset:
                return False
        if withEnd:
            end = max(multimode([frag.interval.end for frag in self.group[index]]))
            if abs(end - interval.end) > self.maxOffset:
                return False
        return True

    def _rewind(self, interval):
        return self._sameLocation(0, interval)

    def _noOverlap(self, interval):
        for frag in self.group:
            if frag[0].interval.chr == interval.chr:
                if interval.start <= frag[0].interval.end and interval.end >= frag[0].interval.start:
                    return False
        return True

    def add(self, fragment):
        if len(self.group) == 0:
            self.group.append([fragment])
            return True
        if self.nfrag == 0:
            if self._rewind(fragment.interval):
                self.group[0].append(fragment)
                self.nfrag = len(self.group)
                self.passes = 1
                if self.pointer == 0:
                    self.passes += 1
                self.pointer = 0
                return True
            else:
                if self._noOverlap(fragment.interval):
                    self.group.append([fragment])
                    self.pointer += 1
                    return True
        else:
            pointer = self.pointer + 1
            if pointer == self.nfrag:
                pointer = 0
            if self._sameLocation(pointer, fragment.interval):
                self.group[pointer].append(fragment)
                self.pointer = pointer
                if self.pointer == self.nfrag - 1:
                    self.passes += 1
                return True
        return False

    def addTerminals(self, firstFragment, lastFragment):
        if len(self.group) == 0:
            return
        self.addLast = False
        self.addFirst = False
        # Add last fragment
        if self.nfrag > 0:
            pointer = self.pointer + 1
            if pointer == self.nfrag:
                pointer = 0
        else:
            pointer = 0
        if lastFragment.interval.strand == '+':
            withStart = True
            withEnd = False
        else:
            withStart = False
            withEnd = True
        if self._sameLocation(pointer, lastFragment.interval, withStart, withEnd):
            self.addLast = True
            self.pointer = pointer
        # Add first fragment
        if firstFragment.interval.strand == '+':
            withStart = False
            withEnd = True
        else:
            withStart = True
            withEnd = False
        if self._sameLocation(-1, firstFragment.interval, withStart, withEnd):
            self.addFirst = True

        if self.addLast:
            self.group[pointer].append(lastFragment)
            if self.nfrag == 0:
                self.nfrag = len(self.group)
                self.passes = 1
        if self.addFirst:
            self.group[-1].append(firstFragment)
            if self.nfrag == 0:
                self.nfrag = len(self.group)
                self.passes = 1
        if self.addLast and self.addFirst:
            if pointer == self.nfrag - 1 and firstFragment.interval.strand == lastFragment.interval.strand:
                if firstFragment.interval.strand == '+':
                    if firstFragment.interval.start <= lastFragment.interval.end + 1:
                        self.passes += 1
                else:
                    if lastFragment.interval.start <= firstFragment.interval.end + 1:
                        self.passes += 1

    def reduce(self):
        if self.nfrag == 0:
            return None
        for frags in self.group:
            chr = frags[0].interval.chr
            start = min(multimode([frag.interval.start for frag in frags]))
            end = max(multimode([frag.interval.end for frag in frags]))
            # 链的方向由链的众数确定
            strand = multimode([frag.interval.strand for frag in frags])[0]
            self.intervals.append(Interval(chr, start, end, strand))
        return self

    def __reduce__(self):
        # 为了使用多线程，所以实现序列化方法
        """
        Enable pickling / object passing.
        """
        state = self.__dict__
        return (self.__class__, (self.maxOffset,), state)

    def to_dict(self):
        """ONLY used for multiprocessing return"""
        return {
            "nfrag": self.nfrag,
            "passes": self.passes,
            "intervals": [
                (iv.chr, iv.start, iv.end, iv.strand)
                for iv in self.intervals
            ]
        }

    @staticmethod
    def from_dict(d):
        fg = FragmentGroup()
        fg.nfrag = d["nfrag"]
        fg.passes = d["passes"]
        fg.intervals = [
            Interval(chr, start, end, strand)
            for chr, start, end, strand in d["intervals"]
        ]
        return fg

# 管理一个完整测序读段的所有片段，包括片段排序、重叠 / 间隙校正、低质量过滤，为后续片段分组提供预处理后的片段列表。
class Read(object):
    ###[todo]先定死后续可以再次更改
    def __init__(self, name=None, length=0, minMapQual=30):
    # def __init__(self, name=None, length=0, minMapQual=args.minMapQual):
        self.name = name
        self.length = int(length)
        self.fragments = []
        self.minMapQual = minMapQual

    def isEmpty(self):
        return len(self.fragments) == 0

    def add(self, fragment):
        self.fragments.append(fragment)

    def ordered(self):
        self.fragments.sort(key=lambda x: x.start)
        return self

    def threadFix(self):
        self.ordered()
        for i in range(len(self.fragments)-1):
            dislen = self.fragments[i+1].start - self.fragments[i].end - 1
            if dislen > 0:
                if self.fragments[i].interval.strand == '+':
                    self.fragments[i].end += dislen
                    self.fragments[i].interval.end += dislen
                    # 在cigar值后面添加M，即匹配操作
                    # [Q]此操作合理吗？
                    self.fragments[i].cigar += '{0}M'.format(dislen)
                else:
                    if self.fragments[i+1].interval.strand == '-':
                        self.fragments[i+1].start -= dislen
                        self.fragments[i+1].interval.end += dislen
                        self.fragments[i+1].cigar += '{0}M'.format(dislen)
                    else:
                        self.fragments[i+1].start -= dislen
                        self.fragments[i+1].interval.start -= dislen
                        self.fragments[i+1].cigar = '{0}M'.format(dislen) + self.fragments[i+1].cigar
            elif dislen < 0:
                dislen = -dislen
                if self.fragments[i].interval.strand == '+':
                    self.fragments[i].end -= dislen
                    # [M]基因组特意的部分，cigar值
                    cigar_re = re.compile(r"(\d+)([MDI])")
                    ces = [[int(m.group(1)), m.group(2)] for m in cigar_re.finditer(self.fragments[i].cigar)]
                    p = len(ces) - 1
                    while dislen > 0 and p >= 0:
                        if ces[p][1] == 'M' or ces[p][1] == 'I':
                            s = min(ces[p][0], dislen)
                            dislen -= s
                            ces[p][0] -= s
                            if ces[p][1] == 'M':
                                self.fragments[i].interval.end -= s
                            if ces[p][0] == 0:
                                p -= 1
                        if p >= 0 and ces[p][1] == 'D':
                            self.fragments[i].interval.end -= ces[p][0]
                            ces[p][0] = 0
                            p -= 1
                    if p < 0:
                        return False
                    # [Q]插入直接变匹配，这是合理的吗
                    if ces[p][1] == 'I':
                        ces[p][1] = 'M'
                        self.fragments[i].interval.end += ces[p][0]
                    self.fragments[i].cigar = ''.join([str(step) + op for step, op in ces[0:p+1]])
                else:
                    if self.fragments[i+1].interval.strand == '-':
                        self.fragments[i+1].start += dislen
                        cigar_re = re.compile(r"(\d+)([MDI])")
                        ces = [[int(m.group(1)), m.group(2)] for m in cigar_re.finditer(self.fragments[i+1].cigar)]
                        p = len(ces) - 1
                        while dislen > 0 and p >= 0:
                            if ces[p][1] == 'M' or ces[p][1] == 'I':
                                s = min(ces[p][0], dislen)
                                dislen -= s
                                ces[p][0] -= s
                                if ces[p][1] == 'M':
                                    self.fragments[i+1].interval.end -= s
                                if ces[p][0] == 0:
                                    p -= 1
                            if p >= 0 and ces[p][1] == 'D':
                                self.fragments[i+1].interval.end -= ces[p][0]
                                ces[p][0] = 0
                                p -= 1
                        if p < 0:
                            return False
                        if ces[p][1] == 'I':
                            ces[p][1] = 'M'
                            self.fragments[i+1].interval.end += ces[p][0]
                        self.fragments[i+1].cigar = ''.join([str(step) + op for step, op in ces[0:p+1]])
                    else:
                        self.fragments[i+1].start += dislen
                        cigar_re = re.compile(r"(\d+)([MDI])")
                        ces = [[int(m.group(1)), m.group(2)] for m in cigar_re.finditer(self.fragments[i+1].cigar)]
                        p = 0
                        while dislen > 0 and p < len(ces):
                            if ces[p][1] == 'M' or ces[p][1] == 'I':
                                s = min(ces[p][0], dislen)
                                dislen -= s
                                ces[p][0] -= s
                                if ces[p][1] == 'M':
                                    self.fragments[i+1].interval.start += s
                                if ces[p][0] == 0:
                                    p += 1
                            if p < len(ces) and ces[p][1] == 'D':
                                self.fragments[i+1].interval.start += ces[p][0]
                                ces[p][0] = 0
                                p += 1
                        if p >= len(ces):
                            return False
                        if ces[p][1] == 'I':
                            ces[p][1] = 'M'
                            self.fragments[i+1].interval.start -= ces[p][0]
                        self.fragments[i+1].cigar = ''.join([str(step) + op for step, op in ces[p:]])
        return True

    def bootstrap(self):
        # 过滤操作，小于3个的都进行过滤
        if len(self.fragments) < 3:
            return None
        res = self.threadFix()
        if res == False:
            return None
        fraggroup = FragmentGroup()
        for i in range(1, len(self.fragments)-1):
            # 若片段比对质量低于阈值，终止循环（不再处理后续片段）
            if self.fragments[i].mapQual < self.minMapQual:
                break
            res = fraggroup.add(self.fragments[i])
            if res == False:
                break
        fraggroup.addTerminals(self.fragments[0], self.fragments[-1])
        fraggroup = fraggroup.reduce()
        return fraggroup


def trans_form_rawalign(read_id, raw, hit_list, mapq):

    read = Read(name=read_id, length=raw, minMapQual=mapq)
    
    for hit in hit_list:
        frag = Fragment(start=hit.q_st, end=hit.q_en,
                        chr=hit.ctg, chrStart=hit.r_st, chrEnd=hit.r_en, strand=hit.strand,
                        mapQual=hit.mapq, cigar=hit.cigar_str)
        read.add(frag)

    return read


def raw_get_ecc(read_id, raw, raw_hit_list, mapq):  # 注意这里添加了minimap_aligner参数，确保能正确引用

    read = trans_form_rawalign(read_id, len(raw), raw_hit_list, mapq)
    fraggroup = read.bootstrap()

    return fraggroup, read


def ccs_hit_to_fraggroup(segments, ccs_hit_list):
    """
    将CCS比对结果(ccs_hit_list)转换为与raw_fraggroup一致的FragmentGroup格式
    """

    fraggroup = FragmentGroup()

    fraggroup.nfrag = len(ccs_hit_list)
    fraggroup.passes = segments.count(";") + 1

    for hit in ccs_hit_list:

        # 创建Fragment对象：参数顺序/含义与raw比对的Fragment完全一致
        frag = Fragment(
            start=hit.q_st,          #  query序列（CCS）上的起始位置
            end=hit.q_en,            #  query序列（CCS）上的结束位置
            chr=hit.ctg,             #  参考基因组contig名称
            chrStart=hit.r_st,       #  参考基因组上的起始位置
            chrEnd=hit.r_en,         #  参考基因组上的结束位置
            strand=hit.strand,       #  链方向（+1为正链，-1为负链，需与raw_hit_list保持一致）
            mapQual=hit.mapq,        #  比对质量（CCS比对质量通常较高，默认100，可按需调整）
            cigar=hit.cigar   #  CIGAR字符串（记录比对细节）
        )

        fraggroup.group.append([frag]) # 将Fragment添加到Read的片段列表中

        fraggroup.intervals.append(Interval(hit.ctg, hit.r_st, hit.r_en, hit.strand))

    return fraggroup
