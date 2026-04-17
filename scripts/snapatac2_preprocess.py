import snapatac2 as snap
import argparse
import os
import sys 
import pandas as pd


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
    "--level_indir", 
    type = str)
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
    "--in_tissue", 
    type = bool,
    default = True)
parser.add_argument(
    "--tissue_barcode",
    type = str,
    default = None)
parser.add_argument(
    "--count", 
    action='store_true')

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
    in_tissue,
    count,
    level_indir,
    n_features,
    max_iter,
    _bin_size,
    tissue_barcode
):
    """
        preprocess bam file and returns fragments and metadata file
    """

    ### out file 
    if in_tissue:
        sampleid_prefix = sampleid 
    else:
        sampleid_prefix = sampleid+ "_all_spots"

    out_fragment = os.path.join(output_dir, sampleid_prefix + "_" + barcode_tag + ".fragments.tsv.gz")
    out_tss = os.path.join(output_dir, sampleid_prefix + "_" + barcode_tag + ".tss_enrichment.pdf")
    out_gene_matrix_h5 = os.path.join(output_dir, sampleid_prefix + "_" + barcode_tag + ".gene_activity.h5ad")

    if _bin_size == "2000":
        out_file_h5 = os.path.join(output_dir, sampleid_prefix + "_" + barcode_tag +".h5ad")
        out_metadata = os.path.join(output_dir, sampleid_prefix + "_" + barcode_tag + ".metadata.txt")
        outpdf = os.path.join(output_dir, sampleid_prefix + "_" + barcode_tag + ".cluster.pdf")
    else:
        out_file_h5 = os.path.join(output_dir, sampleid_prefix + "_" + barcode_tag + "." + _bin_size +".h5ad")
        out_metadata = os.path.join(output_dir, sampleid_prefix + "_" + barcode_tag + "." + _bin_size + ".metadata.txt")
        outpdf = os.path.join(output_dir, sampleid_prefix + "_" + barcode_tag + "." + _bin_size + ".cluster.pdf")

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
    snap.metrics.tsse(data, gene_anno = anno_gff[species])
    data.obs.to_csv(out_metadata.replace("metadata", "all_metadata"), sep="\t",header=True)

    if in_tissue:
        if tissue_barcode is None:
            tag_num = barcode_tag.replace("L", "")
            if barcode_tag == "LB":
                tag_num = "13"
            in_tissue_barcode_file = os.path.join(level_indir, "level_" + tag_num, "barcodes_in_tissue.tsv.gz")
        else:
            in_tissue_barcode_file = tissue_barcode

        print(f'Only processing in-tissue barcode from {in_tissue_barcode_file}', file=sys.stderr)
        print(f'Original barcode: {len(data.obs_names)}', file=sys.stderr)

        cells = pd.read_csv(in_tissue_barcode_file, header=None)
        data = data[data.obs_names.isin(list(cells[0]))]
        print(f'In-tissue barcode: {len(data.obs_names)}', file=sys.stderr)

    # snap.metrics.tsse(data, gene_anno = anno_gff[species])

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
    # if barcode_tag in ["L5", "L6", "L7", "LB", "CB", "SC"]:
    #     data.obs['batch'] = sampleid + "_" + barcode_tag
    #     snap.ex.export_coverage(data, 
    #         groupby='batch',
    #         out_dir=output_dir,
    #         blacklist=blacklist[species],
    #         prefix=sampleid + "_" + barcode_tag + "_bulk",
    #         n_jobs=-1)
    #     ## cluster
    #     snap.ex.export_coverage(data, 
    #         groupby='leiden',
    #         out_dir=output_dir,
    #         blacklist=blacklist[species],
    #         prefix=sampleid + "_" + barcode_tag,
    #         n_jobs=-1)

    # gene_acivity = snap.pp.make_gene_matrix(data, anno_gff[species])
    # gene_acivity.write(out_gene_matrix_h5, compression="gzip")
    
    # data.close()
    print(f'{barcode_tag} has finished !', file=sys.stderr)


def main():
    args = parser.parse_args()
    print(f'is_paired: {args.is_paired}', file=sys.stderr)
    print(f'build: {args.species}', file=sys.stderr)
    print(f'gene_anno: {anno_gff[args.species]}', file=sys.stderr)
    if args.level_indir:
        tags_level =  ["L" + str(i) for i in range(1, 8,1) ]
        tags_level.append("LB")
        print(f'level matrix: {args.level_indir}', file=sys.stderr)
        for tag in tags_level[::-1]:
            print(f'{tag} processing...', file=sys.stderr)
            # for _bin_size in ["500", "2000", "5000"]:
            for _bin_size in ["2000"]:
                bam2fragments(
                    args.input_bam,
                    args.output_dir,
                    args.sampleid,
                    args.is_paired,
                    tag,
                    args.species,
                    args.min_fragments_n,
                    args.max_fragments_n,
                    args.min_tsse,
                    args.in_tissue,
                    args.count,
                    args.level_indir,
                    args.n_features,
                    args.max_iter,
                    _bin_size,
                    args.tissue_barcode
                )
                if args.in_tissue:
                    bam2fragments(
                        args.input_bam,
                        args.output_dir,
                        args.sampleid,
                        args.is_paired,
                        tag,
                        args.species,
                        args.min_fragments_n,
                        args.max_fragments_n,
                        args.min_tsse,
                        False,
                        True,
                        args.level_indir,
                        args.n_features,
                        args.max_iter,
                        _bin_size,
                        args.tissue_barcode
                    )
    else:
        print(f'barcode_tag: {args.barcode_tag}', file=sys.stderr)
        if args.barcode_tag == "SC":
            bin_size_list = ["2000"]
        else:
            bin_size_list = ["500", "2000", "5000"]
        for _bin_size in bin_size_list:
            bam2fragments(
                args.input_bam,
                args.output_dir,
                args.sampleid,
                args.is_paired,
                args.barcode_tag,
                args.species,
                args.min_fragments_n,
                args.max_fragments_n,
                args.min_tsse,
                args.in_tissue,
                args.count,
                None,
                args.n_features,
                args.max_iter,
                _bin_size,
                args.tissue_barcode
            )


if __name__ == '__main__':
    main()