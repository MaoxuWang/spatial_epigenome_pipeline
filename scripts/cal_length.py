#!/usr/bin/env python3

"""
高效率、并行地计算 BAM 文件中每个 CB 标签（细胞）的
read 比对宽度（alignment width）的中位数。

[V2 - 修复] 此版本修复了 'multiprocessing' 的 'Cant pickle local object' 错误。

用法:
    python3 get_read_widths.py \\
        -b /path/to/your/possorted_genome_bam.bam \\
        -t CB \\
        -o /path/to/output/median_lengths.tsv \\
        -c 16
"""

import pysam
import argparse
from collections import defaultdict
import sys
import os
import numpy as np
import multiprocessing
from functools import partial

# 我们在*顶层*定义一个函数来创建嵌套的 defaultdict。
# 这可以被 'pickle' 库序列化，解决了 multiprocessing 的错误。
def nested_int_defaultdict():
    """返回一个 'defaultdict(int)'。"""
    return defaultdict(int)

def find_median_from_counts(lengths_dict):
    """
    内存高效的中位数计算器。
    它通过排序和累加计数来找到中位数，
    而不是在内存中创建一个巨大的列表。
    """
    total_reads = sum(lengths_dict.values())
    if total_reads == 0:
        return 0

    # 将 (length, count) 键值对排序
    sorted_lengths = sorted(lengths_dict.items())
    
    cumulative_count = 0
    
    # 奇数个 reads
    if total_reads % 2 == 1:
        median_index = total_reads // 2
        for length, count in sorted_lengths:
            if cumulative_count + count > median_index:
                return length
            cumulative_count += count
    
    # 偶数个 reads
    else:
        idx1 = (total_reads // 2) - 1
        idx2 = total_reads // 2
        median1 = None
        median2 = None
        
        for length, count in sorted_lengths:
            if median1 is None and cumulative_count + count > idx1:
                median1 = length
            if median2 is None and cumulative_count + count > idx2:
                median2 = length
                # 两个都找到了，可以提前退出
                return (median1 + median2) / 2.0
            cumulative_count += count
    return 0 # 理论上不应到达这里

def process_chunk(chrom_list, bam_file_path, tag):
    """
    (工作进程函数)
    处理一个染色体列表，返回一个计数字典。
    """
    
    # 不再使用 lambda，而是使用我们在顶层定义的 'nested_int_defaultdict'
    local_cell_lengths = defaultdict(nested_int_defaultdict)
    
    try:
        # 每个工作进程必须打开自己的文件句柄
        with pysam.AlignmentFile(bam_file_path, "rb") as bamfile:
            for chrom in chrom_list:
                # fetch 会遍历该染色体上的所有 reads
                for read in bamfile.fetch(chrom):
                    if read.has_tag(tag):
                        barcode = read.get_tag(tag)
                        # read.reference_length 是从 CIGAR 计算的比对宽度
                        length = read.reference_length
                        if length is not None:
                            local_cell_lengths[barcode][length] += 1
    except Exception as e:
        print(f"工作进程出错: {e}", file=sys.stderr)
        return None
        
    return local_cell_lengths

def merge_results(results_list):
    """
    将来自所有工作进程的计数字典合并为一个全局字典。
    """
    print("--- 正在合并所有核心的结果... ---", file=sys.stderr)
    global_cell_lengths = defaultdict(nested_int_defaultdict) # 同样使用新函数
    
    for local_counts in results_list:
        if local_counts is None:
            continue
        for barcode, lengths_dict in local_counts.items():
            for length, count in lengths_dict.items():
                global_cell_lengths[barcode][length] += count
                
    return global_cell_lengths

def calculate_and_write_medians(global_counts, output_file):
    """
    遍历全局字典，计算每个 barcode 的中位数并写入文件。
    """
    print(f"--- 正在计算中位数并写入: {output_file} ---", file=sys.stderr)
    
    with open(output_file, 'w') as f:
        f.write("barcode\ttotal_reads\tmedian_alignment_length\n")
        
        # 遍历所有 barcodes
        for barcode, lengths_dict in global_counts.items():
            total_reads = sum(lengths_dict.values())
            median_length = find_median_from_counts(lengths_dict)
            f.write(f"{barcode}\t{total_reads}\t{median_length}\n")
            
    print("--- 完成! ---", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description="高效率、并行地计算 BAM 文件中每个 CB 标签（细胞）的 read 比对宽度的中位数。"
    )
    parser.add_argument(
        "-b", "--bam_file", 
        type=str, 
        required=True,
        help="Input BAM file (必须已索引, .bai)"
    )
    parser.add_argument(
        "-t", "--tag", 
        type=str, 
        default="CB",
        help="包含细胞 barcode 的 BAM 标签 (例如 'CB')"
    )
    parser.add_argument(
        "-o", "--output_file", 
        type=str, 
        required=True,
        help="输出的 TSV 文件路径 (例如 'median_lengths.tsv')"
    )
    parser.add_argument(
        "-c", "--cores", 
        type=int, 
        default=os.cpu_count(),
        help="用于处理的 CPU 核心数 (默认: 全部可用核心)"
    )
    
    args = parser.parse_args()

    print(f"--- 正在使用 {args.cores} 个核心 ---", file=sys.stderr)
    try:
        with pysam.AlignmentFile(args.bam_file, "rb") as bamfile:
            # 获取所有已比对的染色体
            chromosomes = [ref['SN'] for ref in bamfile.header['SQ'] if bamfile.get_reference_length(ref['SN']) > 0]
    except Exception as e:
        print(f"错误: 无法读取 BAM 文件或其索引: {e}", file=sys.stderr)
        print("请确保 BAM 文件路径正确，并且 .bai 索引文件存在于同一目录。", file=sys.stderr)
        sys.exit(1)
    
    if not chromosomes:
        print("错误：无法从 BAM header 获取染色体列表。", file=sys.stderr)
        sys.exit(1)
        
    # 将染色体列表分成 N 块 (每个核心一块)
    # np.array_split 即使在染色体数 < 核心数时也能正常工作
    chromosome_chunks = np.array_split(chromosomes, args.cores)
    # 转换为 Python 列表
    chromosome_chunks = [chunk.tolist() for chunk in chromosome_chunks]

    print(f"--- 已将 {len(chromosomes)} 个染色体分配给 {len(chromosome_chunks)} 个工作进程 ---", file=sys.stderr)

    # functools.partial 用于“预打包” 'process_chunk' 函数的固定参数
    worker_func = partial(process_chunk, 
                          bam_file_path=args.bam_file, 
                          tag=args.tag)
    
    with multiprocessing.Pool(args.cores) as pool:
        # pool.map 会在所有工作进程上运行 worker_func
        # 并等待所有结果返回
        results_list = pool.map(worker_func, chromosome_chunks)
    
    merged_counts = merge_results(results_list)
    
    calculate_and_write_medians(merged_counts, args.output_file)

if __name__ == "__main__":
    main()