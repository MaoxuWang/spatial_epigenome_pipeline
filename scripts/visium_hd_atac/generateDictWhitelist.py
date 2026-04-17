import pysam
from dscHiCtools import utils 
import os 
import multiprocess
import subprocess
from dscHiCtools.chunkFiles import chunk_bam
from collections import defaultdict, Counter
import gzip
from argparse import ArgumentParser


def getDictChunk(intervals, input_bam):
    """
        extract CR and CB tag from bam file
        Args:
            intervals: slice of bam file that is used for parallel computing
            input_bam: input bam file
        Returns:
            dictionary of CB:CR count
    """
    samfile = pysam.AlignmentFile(input_bam, "rb")
    header = samfile.header
    true_barcode_counter = defaultdict(Counter) 

    for i in intervals:
        for read in samfile.fetch(i[0], i[1], i[2]):
            try:
                true_barcode = read.get_tag('CB')
                UMI = read.get_tag('CB')
            except Exception as e:
                continue
            true_barcode_counter[true_barcode][raw_barcode] += 1

    samfile.close()
    return true_barcode_counter


@utils.log_info
def process_whitelist(args):
    """
        add cell barcode tags to bam file
        input:
            bam file, the header of read name is barcode 
        output:
            bam file, with "CB:" tag added
    """
    inputBam = pysam.AlignmentFile(args.input_bam, "rb")
    idx_file = args.input_bam + ".bai"
    if not os.path.exists(idx_file):
        utils.eprint("[Whitelist::] Indexing input bam file...")
        Args_m = [f'{args.samtools} index -@ {args.threads} {args.input_bam}']
        subprocess.check_call(Args_m, shell=True)
    intervals = chunk_bam(inputBam, args.threads)

    inputBam.close()
    p = multiprocess.Pool(args.threads)

    utils.eprint("[Whitelist::] create barcode whitelist...")
    utils.eprint(f"[Whitelist::] Multi-cores: {args.threads}")
    parallel_results = []

    for interval in intervals.values():
        p.apply_async(
            getDictChunk,
            args=(
                interval,
                args.input_bam,
            ),
            error_callback=utils.print_error,
            callback=parallel_results.append,
        )
    p.close()
    p.join()
    utils.eprint("[Whitelist::] Filtering done.")
    utils.eprint("[Whitelist::] Merging results...")
    
    merged_counter = mergeChunks(parallel_results)
    generateWhitelist(merged_counter, args.outdir)

def parse_args(parser):
    parser.add_argument(
        '-i', '--input_bam', 
        required = True
    )
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
