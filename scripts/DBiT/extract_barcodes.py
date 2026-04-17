import argparse
import pybktree
import Levenshtein
from dscHiCtools import utils 
import os 
import sys 
import time 
from collections import Counter
import gzip 
import multiprocess
import subprocess
from itertools import islice
from dscHiCtools.write import (
    write_unmapped, 
    write_mapped
)
from dscHiCtools.chunkFiles import (
    get_n_lines, 
    chunk_reads
)
from argparse import ArgumentParser
from pathlib import Path

def parse_tuple(s):
    try:
        # 去掉可能存在的括号，并按逗号分割
        # 例如 "(98, 108)" -> "98, 108" -> ["98", " 108"]
        items = s.strip('()').split(',')
        # 转换为整数元组
        return tuple(int(x.strip()) for x in items)
    except:
        raise argparse.ArgumentTypeError(f"格式错误: {s}。请使用类似 '(start, end)' 的格式。")
                                         
## mapping:     
def makeBKtree(barcode_f):
    """
        make the BKtree for reference barcode file;
        This is for fast matching
        Returns:
            a BKtree
    """
    barcode_tree = pybktree.BKTree(Levenshtein.hamming, barcode_f)
    utils.eprint("[mapBarcode::] Generated barcode tree from reference barcode file")
    return(barcode_tree)


def extractBarcode(
    seq,
    cb_array,
    umi_array,
):
    """
    """
    cell_barcodes = []
    UMI_barcodes = []

    for pos in cb_array:
        cell_barcodes.append(seq[pos[0]:pos[1]])
    for pos in umi_array:
        UMI_barcodes.append(seq[pos[0]:pos[1]])
    return "".join(cell_barcodes), UMI_barcodes[0]


def writeFile(
    outFile,
    cellBarcode,
    UMI,
    readNames,
    seqs,
    chars,
    qualities,
    i,
    paired
):
    # readName = readNames[i].strip().replace(" ", "_") + ":" + cellBarcode + ":" + UMI + "\n"  
    # readName = readNames[i].strip() + ":" + cellBarcode + ":" + UMI + "\n"  
    readName = "@" + cellBarcode + ":" + UMI + ":" + readNames[i][1:]  

    outFile.write(readName.encode())

    if paired :
        outFile.write(seqs[i][79:].encode())
        outFile.write(chars[i].encode())
        outFile.write(qualities[i][79:].encode())
    else:
        outFile.write(seqs[i].encode())
        outFile.write(chars[i].encode())
        outFile.write(qualities[i].encode())


def autoDetectChemistry(
    read2,
    barcode_tree,
    cb_array,
    umi_array,
    min_mismatch,
    n = 40000):
    """
        auto detect whether to reverse the read2 cell barcode
        Args:
            read2: R2 file
            barcode_tree: BKtree for reference barcode file
            min_mismatch: min mismatch bases allowed for mapping cell barcode to reference
        Returns: 
            reverse: whether to reverse barcode (associated with chemistry)
    """
    utils.eprint("[main::] Auto detecting orientation...")
    i = 0
    fd_result =0
    rs_result = 0

    with gzip.open(read2, 'rt') as textfile1:
        totallines = islice(
            textfile1, 0 * 4, n * 4
        )
        while True:
            try:    
                readNames = next(totallines)
                seq = next(totallines).strip()
                _char = next(totallines)
                qualities = next(totallines)

            except StopIteration as e:
                break

            cell_barcode_string, _ = extractBarcode(seq, cb_array, umi_array)
            cell_barcode_string_reverse = utils.reverseRead(cell_barcode_string)
            if utils.alignReference(cell_barcode_string, barcode_tree, min_mismatch) != 0:
                fd_result += 1
            if utils.alignReference(cell_barcode_string_reverse, barcode_tree, min_mismatch) != 0:
                rs_result += 1
            i += 1

    utils.eprint(f'[dscHiCtools::] Forward: {fd_result} hits of {i} reads ({round(fd_result / i, 4) * 100}% rate)')
    utils.eprint(f'[dscHiCtools::] Reverse: {rs_result} hits of {i} reads ({round(rs_result / i, 4) * 100}% rate)')
    if fd_result > rs_result:
        utils.eprint("[dscHiCtools::] orientation sets to be forward")
        return False
    else:
        utils.eprint("[dscHiCtools::] orientation sets to be reverse")
        return True 


def mapCellBarcode(
    read1,
    read2,
    cb_array,
    umi_array,
    indexes,
    out_r1,
    out_r2,
    reverse,
    barcode_tree,
    min_mismatch,
    paired
):
    """
        Add cell barcodes to reads through indexes of input file (for multi-thread running)
        Args:
            input: input R1.fastq R2.fastq file containing barcode reads
            out_r1: output ; read name added by cell barcode
            out_r2
            barcode_tree: BKtree for reference barcode file
            min_mismatch: min mismatch bases allowed for mapping cell barcode to reference
        Returns:
            mapped.mtx: mapped cell barcode - readsN
            unmapped.mtx: unmapped cell barcode - readsN - sequence
            barcoded reads
    """
    ## count the processed reads 
    n = 1
    t = time.time()
    mapped = Counter()
    unmapped = Counter()

    ## must output to a .gz file
    out_r1_file = gzip.open(out_r1, "wb")
    out_r2_file = gzip.open(out_r2, "wb")

    with gzip.open(read1, "rt") as textfile1, gzip.open(read2, "rt") as textfile2:
        totallines = islice(
            zip(textfile1, textfile2), indexes[0] * 4, indexes[1] * 4
        )
        while True:
            try:
                readNames = next(totallines)
                seqs = next(totallines)
                _char = next(totallines)
                qualities = next(totallines)

                # Progress info
                if n % 10000000 == 0:
                    utils.eprint(
                        "[mapBarcode::] Processed 10,000,000 reads in {}. Total "
                        "reads: {:,} in child {}".format(
                            utils.secondsToText(time.time() - t), n, os.getpid()
                        )
                    )
                    sys.stdout.flush()
                    t = time.time()
                n += 1
                ## barcode is on read2
                read = seqs[1].strip()

                cell_barcode_string_raw, UMI_raw = extractBarcode(read, cb_array, umi_array)
                if reverse:
                    cell_barcode_string = utils.reverseRead(cell_barcode_string_raw)
                    UMI_string = utils.reverseRead(UMI_raw)
                else:
                    cell_barcode_string = cell_barcode_string_raw
                    UMI_string = UMI_raw
                
                if len(cell_barcode_string) != 16:
                    utils.eprint("[mapBarcode::] wrong cell barcode length found!")
                    continue

                ### unmapped
                cell_barcode = utils.alignReference(cell_barcode_string, barcode_tree, min_mismatch)
                if( not cell_barcode):
                    unmapped[cell_barcode_string] = unmapped.get(cell_barcode_string, 0) + 1 
                    continue
                ### mapped
                else:
                    mapped[cell_barcode] = mapped.get(cell_barcode, 0) + 1 

                    writeFile(out_r1_file, 
                        cell_barcode,
                        UMI_string,
                        readNames,
                        seqs,
                        _char,
                        qualities,
                        0,
                        False)
                    writeFile(out_r2_file, 
                        cell_barcode,
                        UMI_string,
                        readNames,
                        seqs,
                        _char,
                        qualities,
                        1,
                        paired)
                   
            except StopIteration as e:
                break
            except Exception as e:
                utils.eprint(f"[Barcode::] {str(e)}")
                continue
    
    out_r1_file.close()
    out_r2_file.close()
    utils.eprint(
        "[mapBarcode::] Mapping done for process {}. Processed {:,} reads".format(os.getpid(), n - 1)
    )
    sys.stdout.flush()
    return mapped, unmapped

def checkFile(file_list):
    for _file in file_list:
        if not os.path.exists(_file) or os.path.getsize(_file) == 0:
            raise FileNotFoundError(f"File {_file} is empty ! Please re-run.")

@utils.log_info
def extract_barcodes(args):
    """
        add cell barcode to read name of R1, R3
        Args:
            read1: read1 file from 10X data (.fastq.gz)
            read2: read2 file from 10X data (.fastq.gz)
            outdir: directory for output
            threads: >= 1
            sampleName: output prefix that specifys sample
            ref_barcode: reference cell barcode file (.txt) 
            min_mismatch: minimal mismatch allowed for cell barcode mapping
            n_lines: file line count
        Returns:
            barcoded.fastq.gz (using RNA barcode)
            QC/ mapped.tsv.gz: mapped cell barcode - readsN
            QC/ unmapped.tsv.gz: unmapped cell barcode - readsN - sequence
    """
    if args.n_lines is None or args.n_lines <= 0 :
        args.n_lines = get_n_lines(args.read2)
    n_reads = int(args.n_lines / 4)
    utils.eprint("[mapBarcode::] Start adding cell barcodes...")
    utils.eprint(f"[mapBarcode::] Total input fastq files have {args.n_lines} lines.\nProcessing {n_reads} reads")

    with open(os.path.join(args.outdir, args.sampleName + ".total_readsN.txt"), 'w') as readN_f:
        readN_f.write(str(n_reads))
    
    DNA_barcodes = utils.readBarcode(args.ref_barcode)
    barcode_tree = makeBKtree(DNA_barcodes)
    if not os.path.isfile(args.read1) :
        raise Exception(f"{args.read1} not found !")

    if not os.path.isfile(args.read2) :
        raise Exception(f"{args.read2} not found !")

    ## matrix out to QC subdirectory
    outfolder = os.path.join(os.path.dirname(args.outdir), "QC")

    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    ## detect strand orientation
    reverse = autoDetectChemistry(
            args.read2,
            barcode_tree,
            args.cb_array,
            args.umi_array,
            args.min_mismatch
        )

    if args.threads <= 1:
        out_r1 = os.path.join(args.outdir, args.sampleName + ".R1.barcoded.fastq.gz")
        out_r2 = os.path.join(args.outdir, args.sampleName + ".R2.barcoded.fastq.gz")

        mappedMtx, unmappedMtx = mapCellBarcode(
            read1=args.read1,
            read2=args.read2,
            cb_array=args.cb_array,
            umi_array=args.umi_array,
            indexes=[0, n_reads],
            out_r1=out_r1,
            out_r2=out_r2,
            reverse=reverse,
            barcode_tree=barcode_tree,
            min_mismatch=args.min_mismatch,
            paired=args.paired
        )
        utils.eprint("[mapBarcode::] Cell barcodes are added.")
    else:
        utils.eprint(f"[mapBarcode::] Adding cell barcodes is running with {args.threads} cores.")
        p = multiprocess.Pool(processes=args.threads)
        chunk_indexes = chunk_reads(n_reads, args.threads)

        r1_single_files = []
        r2_single_files = []
        parallel_results = []

        ## multi-processing
        for i, indexes in enumerate(chunk_indexes):
            out_r1 = os.path.join(args.outdir, args.sampleName + "_" + str(i) + ".R1.barcoded.fastq.gz")
            out_r2 = os.path.join(args.outdir, args.sampleName + "_" + str(i) + ".R2.barcoded.fastq.gz")
            r1_single_files.append(out_r1)
            r2_single_files.append(out_r2)

            p.apply_async(
                mapCellBarcode,
                args=(
                    args.read1,
                    args.read2,
                    args.cb_array,
                    args.umi_array,
                    indexes,
                    out_r1,
                    out_r2,
                    reverse,
                    barcode_tree,
                    args.min_mismatch,
                    args.paired
                ),
                error_callback=utils.print_error,
                callback=parallel_results.append,
            )

        p.close()
        p.join()
        utils.eprint("[mapBarcode::] Merging results...")
        mappedMtx, unmappedMtx = utils.merge_results(parallel_results)

        ## merge fastq.gz subprocess
        out_final_r1 = os.path.join(args.outdir, args.sampleName + ".R1.barcoded.fastq.gz")
        out_final_r2 = os.path.join(args.outdir, args.sampleName + ".R2.barcoded.fastq.gz")

        checkFile(r1_single_files)
        checkFile(r2_single_files)

        Args_m = [f'cat {" ".join(r1_single_files)} > {out_final_r1}']
        subprocess.check_call(Args_m, shell=True)

        Args_m = [f'cat {" ".join(r2_single_files)} > {out_final_r2}']
        subprocess.check_call(Args_m, shell=True)

        ## remove temporary files
        for my_file in r1_single_files:
            if os.path.exists(my_file):
                os.remove(my_file)
        for my_file in r2_single_files:
            if os.path.exists(my_file):
                os.remove(my_file)
        utils.eprint("[mapBarcode::] Cell barcodes are added and mapped!")

    ## save to stat
    unmapped_file = args.sampleName + ".Top500_unmapped.tsv.gz"
    mapped_file = args.sampleName + ".mapped.tsv.gz"

    ### unmapped barcode
    num_unmapped = write_unmapped(
        merged_no_match=unmappedMtx,
        outfolder=outfolder,
        filename=unmapped_file,
    )

    ### mapped barcode
    num_mapped = write_mapped(mappedMtx, outfolder, mapped_file)

    statistic_file = os.path.join(outfolder, args.sampleName + ".barcodes.stat.tsv")
    with open(statistic_file, 'wt') as statistic_f:
        statistic_f.write(f"mapped: {str(num_mapped)}\n")
        statistic_f.write(f"unmapped: {str(num_unmapped)}\n")


def parse_args(parser):
    parser.add_argument('--read1', required = True)
    parser.add_argument('--read2', required = True)
    parser.add_argument('-o', '--outdir', type=str, default="output directory")
    parser.add_argument('--sampleName', type=str, default="out", help='sample name')
    parser.add_argument('-t', '--threads', type=int, default=10)
    parser.add_argument(
                        "--n_lines",
                        type=int,
                        default=-1,
                        help="Total lines of input fastq file to be processed. Default: all (-1)")
    parser.add_argument(
                        "--ref_barcode",
                        type=str,
                        default=str(Path(__file__).resolve().parents[2] / "resources" / "DBiT" / "DBiT.whitelist.txt"),
                        help="The DBiT reference barcode file")
    parser.add_argument('--cb_array', nargs='+', 
        type = parse_tuple, 
        help = "list of position of cell barcodes in Read2",
        default=[(32, 40), (70, 78)]
    )
    parser.add_argument('--umi_array', nargs='+', 
        type = parse_tuple, 
        help = "list of position of UMI in Read2",
        default=[(22, 32)]
    )
    parser.add_argument(
        "--paired",
        action='store_true')
    parser.add_argument(
                        "--min_mismatch",
                        type=int,
                        default=1,
                        help="The minimum mismatch between fastq barcode and ref barcode")



    return parser.parse_args()


def main():
    """
        extract barcode & UMI from R2
    """
    parser = ArgumentParser(description=main.__doc__)
    args = parse_args(parser)

    extract_barcodes(args)

if __name__ == "__main__":
    main()
