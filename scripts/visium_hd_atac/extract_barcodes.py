from dscHiCtools import utils 
import os 
import multiprocess
import subprocess
from collections import defaultdict, Counter
import gzip
from argparse import ArgumentParser
from dscHiCtools.chunkFiles import (
    get_n_lines, 
    chunk_reads
)
from itertools import islice
import time
import sys
import re


def extract_tags(read):
    cb_value = -1  # default no result
    ub_value = None
    
    # search CB:Z
    for field in read:
        if re.match(r'^CB:Z:', field):
            parts = field.split('Z:')
            if len(parts) > 1:
                cb_value = parts[-1].split()[0]
                break 
    
    # search UB:Z
    for field in read:
        if re.match(r'^UB:Z:', field):
            parts = field.split('Z:')
            if len(parts) > 1:
                ub_value = parts[-1].split()[0]
                break
    
    if ub_value is None:
        for field in read:
            if re.match(r'^UR:Z:', field):
                parts = field.split('Z:')
                if len(parts) > 1:
                    ub_value = parts[-1].split()[0]
                    break


    return cb_value, ub_value


def extract_read1(read):
    """
    从一个已分割的 SAM/BAM read 记录中提取 1R:Z: 和 1Y:Z: 标签的值。

    Args:
        read (list): 一个字符串列表，代表 SAM/BAM 文件中的一行记录，已按制表符分割。

    Returns:
        tuple: 一个包含两个元素的元组 (r1_value, y1_value)。
               如果某个标签未找到，其对应的值将为 None。
    """
    r1_value = None
    y1_value = None

    # 遍历 read 记录中的每一个字段（标签）
    for field in read:
        # 匹配 "1R:Z:" 标签
        if r1_value is None and re.match(r'^1R:Z:', field):
            try:
                # 分割字符串以提取标签值
                r1_value = field.split('Z:')[1].split()[0]
            except IndexError:
                # 以防万一 "Z:" 后面没有内容
                r1_value = ""
            continue # 找到后继续寻找下一个，避免重复匹配

        # 匹配 "1Y:Z:" 标签
        if y1_value is None and re.match(r'^1Y:Z:', field):
            try:
                # 分割字符串以提取标签值
                y1_value = field.split('Z:')[1].split()[0]
            except IndexError:
                y1_value = ""

        # 优化：如果两个标签都已找到，则提前结束循环
        if r1_value is not None and y1_value is not None:
            break
            
    return r1_value, y1_value


def convertChunk(indexes, input_bam, out_fq, paired):
    """
        extract CR and CB tag from bam file
        Args:
            indexes: slice of bam file that is used for parallel computing
            input_bam: input sam.gz file
        Returns:
            dictionary of CB:CR count
    """
    n = 1
    t = time.time()
    mapped = 0

    textfile1 = gzip.open(input_bam, "rt")
    out_fq_file = gzip.open(out_fq, "wb")
    if paired:
        out_fq_R1 = gzip.open(out_fq.replace("R2", "R1"), "wb")

    totallines = islice(
        textfile1, indexes[0] * 1, indexes[1] * 1
    )
    while True:
        try:
            read = next(totallines).rstrip().split("\t")
            n += 1
            true_barcode, UMI = extract_tags(read)
            if true_barcode == -1:
                continue
            else:
                mapped += 1

            new_name = "@" + true_barcode.replace("-1", "") + ":" + UMI + ":" + read[0]
            if n % 1000000 == 0:
                utils.eprint(
                    "[Barcode::] Processed 1,000,000 reads in {}. Total "
                    "reads: {:,} in child {}".format(
                        utils.secondsToText(time.time() - t), n, os.getpid()
                    )
                )
                sys.stdout.flush()
                t = time.time()
            sequence = read[9]
            qualities = read[10]
            read_info = new_name + "\n" + sequence + "\n+\n" + qualities + "\n"

            if paired:
                r1_seq, r1_qual = extract_read1(read)
                r1_seq = r1_seq[43:]
                r1_qual = r1_qual[43:]
                r1_info = new_name + "\n" + r1_seq + "\n+\n" + r1_qual + "\n"
                out_fq_R1.write(r1_info.encode())

            out_fq_file.write(read_info.encode())
        except StopIteration as e:
            break
        except Exception as e:
            print(read)
            utils.eprint(f"[Barcode::] {str(e)}")
            continue
    utils.eprint(
        "[Barcode::] Done for process {}. Processed {:,} reads, {} valid barcodes".format(os.getpid(), n - 1, mapped)
    )
    sys.stdout.flush()

    return mapped

@utils.log_info
def process_whitelist(args):
    """
        add cell barcode tags to bam file
        input:
            bam file, the header of read name is barcode 
        output:
            bam file, with "CB:" tag added
    """
    if args.n_lines is None or args.n_lines <= 0 :
        args.n_lines = get_n_lines(args.input_bam, _format = "sam")
    ## input: sam.gz ; one line -> one read
    n_reads = int(args.n_lines / 1)
    utils.eprint("[Barcode::] Start adding cell barcodes...")
    utils.eprint(f"[Barcode::] Total input fastq files have {n_reads} reads")

    with open(os.path.join(args.outdir, args.sample + ".total_readsN.txt"), 'w') as readN_f:
        readN_f.write(str(n_reads))

    chunk_indexes = chunk_reads(n_reads, args.threads)

    p = multiprocess.Pool(args.threads)

    utils.eprint("[Barcode::] create barcoded fastq...")
    utils.eprint(f"[Barcode::] Multi-cores: {args.threads}")
    utils.eprint(f"[Barcode::] Chunk indexes: {chunk_indexes}")

    r2_single_files = []
    mapped_result = []

    if args.paired:
        r1_single_files = []

    for i, indexes in enumerate(chunk_indexes):
        # os.makedirs(os.path.join(args.outdir, "tmp"), exist_ok=True)
        # out_fq = os.path.join(args.outdir, "tmp", f"{args.sample}_{i}.barcoded.fastq.gz")
        out_fq = os.path.join(args.outdir, f"{args.sample}_{i}.R2_tmp.fastq.gz")
        r2_single_files.append(out_fq)
        if args.paired:
            out_fq_R1 = out_fq.replace("R2", "R1")
            r1_single_files.append(out_fq_R1)

        p.apply_async(
            convertChunk,
            args=(
                indexes,
                args.input_bam,
                out_fq,
                args.paired
            ),
            error_callback=utils.print_error,
            callback=mapped_result.append,
        )
    p.close()
    p.join()
    utils.eprint("[Barcode::] Merging results...")

    out_final_r2 = os.path.join(args.outdir, args.sample + ".R2.barcoded.fastq.gz")
    for _file in r2_single_files:
        if not os.path.exists(_file) or os.path.getsize(_file) == 0:
            raise FileNotFoundError(f"File {_file} is empty ! Please re-run.")
    Args_m = [f'cat {" ".join(r2_single_files)} > {out_final_r2}']
    subprocess.check_call(Args_m, shell=True)

    ## remove temporary files
    for my_file in r2_single_files:
        if os.path.exists(my_file):
            os.remove(my_file)

    if args.paired:
        out_final_r1 = os.path.join(args.outdir, args.sample + ".R1.barcoded.fastq.gz")
        for _file in r1_single_files:
            if not os.path.exists(_file) or os.path.getsize(_file) == 0:
                raise FileNotFoundError(f"File {_file} is empty ! Please re-run.")
        Args_m = [f'cat {" ".join(r1_single_files)} > {out_final_r1}']
        subprocess.check_call(Args_m, shell=True)

        for my_file in r1_single_files:
            if os.path.exists(my_file):
                os.remove(my_file)

    with open(os.path.join(args.outdir, args.sample + ".barcode_stat.txt"), 'w') as readN_f:
        readN_f.write(f'Total:\t{str(n_reads)}\nValid:\t{str(sum(mapped_result))}')

def parse_args(parser):
    parser.add_argument(
        '-i', '--input_bam', 
        required = True
    )
    parser.add_argument(
        "--n_lines",
        type=int,
        default=-1,
        help="Total lines of input fastq file to be processed. Default: all (-1)")
    parser.add_argument(
        '--samtools', 
        type=str, 
        default=os.environ.get("SAMTOOLS", "samtools")
    )
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=10
    )
    parser.add_argument(
        '--sample',
        type=str,
        required = True
    )
    parser.add_argument(
        "--paired",
        action='store_true',
        help="Enable paired-end mode. If set, read1 sequences will be extracted by 1R 1Y tag and trim first 43bp."
    )
    parser.add_argument(
        "-o", "--outdir",
        type=str,
        required = True,
        help="output directory")
    
    return parser.parse_args()

def main():
    """
        generate whitelist from spaceranger processed bam file
    """
    parser = ArgumentParser(description=main.__doc__)
    args = parse_args(parser)
    process_whitelist(args)

if __name__ == "__main__":
    main()
