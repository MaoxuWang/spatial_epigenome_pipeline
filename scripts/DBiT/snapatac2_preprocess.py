import snapatac2 as snap
import argparse
import os
import sys 
import pandas as pd

"""
    This script is used to generate peak and tile matrix used in spatial ATAC-seq analysis.
    Metadata is generated based on the tissue region and spatial barcodes.
    Only Barcodes of ROI are considered
"""
extra_bin = os.environ.get("SPATIAL_EPIGENOME_EXTRA_BIN")
if extra_bin:
    os.environ["PATH"] = os.environ["PATH"] + os.pathsep + extra_bin

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", "--input_bam",
    required=True,
    type = str)
parser.add_argument(
    "-o", "--output_dir",
    required=True,
    type = str)
parser.add_argument(
    "--sampleid",
    required=True,
    type = str)
parser.add_argument(
    "--is_paired",
    action='store_true')
parser.add_argument(
    "--barcode_tag", 
    type = str,
    default = "CB")
parser.add_argument(
    "--species", 
    type = str)
parser.add_argument(
    "--min_fragments_n", 
    default = 1000,
    type = int)
parser.add_argument(
    "--max_iter", 
    default = 2,
    type = int)
parser.add_argument(
    "--n_features", 
    default = 10000,
    type = int)
parser.add_argument(
    "--max_fragments_n", 
    default = 10000000000,
    type = int)
parser.add_argument(
    "--min_tsse", 
    default = 1,
    type = int)
parser.add_argument(
    "--tissue_barcode",
    type = str,
    default = None)
parser.add_argument(
    "--count", 
    help = "only to generate metadata file",
    action='store_true')
## 500 5000 2000
parser.add_argument(
    "--bin_size_list", 
    nargs = '+',
    default = ["5000"],
    type = str,
    help = "bin size list for tile matrix")


genome = {'mm10': snap.genome.mm10, 'hg38': snap.genome.hg38}
anno_gff = {
    "mm10": os.environ.get("MM10_GFF3", "/path/to/gencode.mm10.annotation.gff3"),
    "hg38": os.environ.get("HG38_GFF3", "/path/to/gencode.hg38.annotation.gff3"),
}
blacklist = {
    "mm10": os.environ.get("MM10_BLACKLIST_BED", "/path/to/ENCODE.mm10.blacklist.bed"),
    "hg38": os.environ.get("HG38_BLACKLIST_BED", "/path/to/ENCODE.hg38.blacklist.bed"),
}


def bam2fragments(
    input_bam,
    output_dir,
    sampleid,
    is_paired,
    barcode_tag,
    species,
    min_fragments_n,
    max_fragments_n,
    min_tsse,
    count,
    n_features,
    max_iter,
    _bin_size,
    tissue_barcode
):
    """
        preprocess bam file and returns fragments and metadata file
    """
    ### out file 
    out_fragment = os.path.join(output_dir, sampleid + "_" + barcode_tag + ".fragments.tsv.gz")
    out_tss = os.path.join(output_dir, sampleid + "_" + barcode_tag + ".tss_enrichment.pdf")
    out_gene_matrix_h5 = os.path.join(output_dir, sampleid + "_" + barcode_tag + ".gene_activity.h5ad")

    # default: 5kb
    if _bin_size == "5000":
        out_file_h5 = os.path.join(output_dir, sampleid + "_" + barcode_tag +".h5ad")
        out_metadata = os.path.join(output_dir, sampleid + "_" + barcode_tag + ".metadata.txt")
        outpdf = os.path.join(output_dir, sampleid + "_" + barcode_tag + ".cluster.pdf")
    else:
        out_file_h5 = os.path.join(output_dir, sampleid + "_" + barcode_tag + "." + _bin_size +".h5ad")
        out_metadata = os.path.join(output_dir, sampleid + "_" + barcode_tag + "." + _bin_size + ".metadata.txt")
        outpdf = os.path.join(output_dir, sampleid + "_" + barcode_tag + "." + _bin_size + ".cluster.pdf")

    snap.pp.make_fragment_file(
        input_bam,
        out_fragment,
        is_paired = is_paired,
        barcode_tag=barcode_tag)
    data = snap.pp.import_data(
            out_fragment,
            chrom_sizes=genome[species],
            sorted_by_barcode=False,
            min_num_fragments=0,
            n_jobs=-1
        )

    print(f'Only processing in-tissue barcode from {tissue_barcode}', file=sys.stderr)
    print(f'Original barcode: {len(data.obs_names)}', file=sys.stderr)

    cells = pd.read_csv(tissue_barcode, header=None)
    data = data[data.obs_names.isin(list(cells[0]))]
    print(f'In-tissue barcode: {len(data.obs_names)}', file=sys.stderr)

    snap.metrics.tsse(data, gene_anno = anno_gff[species])

    snap.pl.tsse(data,
        interactive=False, 
        out_file=out_tss, 
        min_fragment=min_fragments_n)
    if count:
        data.obs.to_csv(out_metadata, sep="\t",header=True)
        return
    ## clustering no filtering
    snap.pp.filter_cells(data, 
        min_counts=min_fragments_n,
        min_tsse=min_tsse,
        max_counts=max_fragments_n)

    snap.pp.add_tile_matrix(data, bin_size = int(_bin_size))
    snap.pp.select_features(data, 
        blacklist=blacklist[args.species],
        n_features=n_features,
        filter_lower_quantile = 0.01,
        filter_upper_quantile = 0.01, 
        max_iter = max_iter)

    snap.tl.spectral(data)
    snap.tl.umap(data)

    snap.pp.knn(data)
    snap.tl.leiden(data)
    snap.pl.umap(data, color='leiden', interactive=False, out_file=outpdf, height=500)
    data.obs.to_csv(out_metadata, sep="\t",header=True)

    data.write(out_file_h5, compression="gzip")
    print(f'{barcode_tag} has finished !', file=sys.stderr)
    return data


def getPeakMatrix(
    atac_st,
    species,
    output_dir,
    sampleid,
    groupby = 'leiden'):

    snap.tl.macs3(
        atac_st, 
        groupby = groupby,
        qvalue=0.01,
        nolambda=True,
        shift=-100,
        extsize=200,
        blacklist=blacklist[args.species]
        )


    peaks = snap.tl.merge_peaks(atac_st.uns['macs3'], genome[species], half_width=250)
    peak_mat = snap.pp.make_peak_matrix(atac_st, use_rep=peaks['Peaks'])
    mat = pd.DataFrame(peak_mat.X.todense(), index=peak_mat.obs_names, columns=peak_mat.var_names)
    out_peak_matrix_h5 = os.path.join(output_dir, sampleid + ".peak_matrix.h5")
    
    mat.to_hdf(out_peak_matrix_h5, key="mat")


def main():
    args = parser.parse_args()
    print(f'is_paired: {args.is_paired}', file=sys.stderr)
    print(f'build: {args.species}', file=sys.stderr)
    print(f'gene_anno: {anno_gff[args.species]}', file=sys.stderr)
  
    print(f'barcode_tag: {args.barcode_tag}', file=sys.stderr)
    for _bin_size in args.bin_size_list:
        atac_st = bam2fragments(
            args.input_bam,
            args.output_dir,
            args.sampleid,
            args.is_paired,
            args.barcode_tag,
            args.species,
            args.min_fragments_n,
            args.max_fragments_n,
            args.min_tsse,
            args.count,
            args.n_features,
            args.max_iter,
            _bin_size,
            args.tissue_barcode
        )
        if not args.count:
            getPeakMatrix(
                atac_st,
                args.species,
                args.output_dir,
                args.sampleid)

if __name__ == '__main__':
    main()