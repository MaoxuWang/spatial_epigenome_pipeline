import pysam
from dscHiCtools import utils 
import os 
import multiprocess
import subprocess
from dscHiCtools.chunkFiles import chunk_bam
from collections import defaultdict, Counter
import gzip
from argparse import ArgumentParser


def convertChunk(intervals, input_bam, out_fq):
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

    out_fq_file = gzip.open(out_fq, "wb")

    for i in intervals:
        for read in samfile.fetch(i[0], i[1], i[2], until_eof=True):
            try:
                true_barcode = read.get_tag('CB')
                UMI = read.get_tag('UB')

            except Exception as e:
                continue
            new_name = "@" + true_barcode.replace("-1", "") + ":" + UMI + ":" + read.query_name
            # ## R1：

            ## R2：
            sequence = read.get_forward_sequence()
            qualities = "".join(chr(q + 33)  for q in read.get_forward_qualities())

            read_info = new_name + "\n" + sequence + "\n+\n" + qualities + "\n"
            out_fq_file.write(read_info.encode())

    samfile.close()
    out_fq_file.close()

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
        utils.eprint("[Barcode::] Indexing input bam file...")
        Args_m = [f'{args.samtools} index -@ {args.threads} {args.input_bam}']
        subprocess.check_call(Args_m, shell=True)
    intervals = chunk_bam(inputBam, args.threads)

    inputBam.close()
    p = multiprocess.Pool(args.threads)

    utils.eprint("[Barcode::] create barcode whitelist...")
    utils.eprint(f"[Barcode::] Multi-cores: {args.threads}")
    r2_single_files = []

    for i, interval in enumerate(intervals.values()):
        out_fq = os.path.join(args.outdir, args.sample + "_" + str(i) + ".barcoded.fastq.gz")
        r2_single_files.append(out_fq)

        p.apply_async(
            convertChunk,
            args=(
                interval,
                args.input_bam,
                out_fq,
            ),
            error_callback=utils.print_error,
        )
    p.close()
    p.join()
    utils.eprint("[Barcode::] Merging results...")
    
    out_final_r2 = os.path.join(args.outdir, args.sample + ".R2.barcoded.fastq.gz")

    Args_m = [f'cat {" ".join(r2_single_files)} > {out_final_r2}']
    subprocess.check_call(Args_m, shell=True)

    ## remove temporary files
    for my_file in r2_single_files:
        if os.path.exists(my_file):
            os.remove(my_file)

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
