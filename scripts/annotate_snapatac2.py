import snapatac2 as snap
import os
import warnings
warnings.filterwarnings("ignore")
import anndata as ad
import pandas as pd
import scanpy as sc
import scvi
import argparse
import sys 
import pickle

scvi.settings.seed = 0
extra_bin = os.environ.get("SPATIAL_EPIGENOME_EXTRA_BIN")
if extra_bin:
    os.environ["PATH"] = os.environ["PATH"] + os.pathsep + extra_bin

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", "--input_h5ad",
    required=True,
    type = str)
parser.add_argument(
    "-r", "--reference_h5ad",
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
    "--species", 
    default = "mm10",
    type = str)
parser.add_argument(
    "--min_fragments_n", 
    default = 1000,
    type = int)
parser.add_argument(
    "--max_fragments_n", 
    default = 100000,
    type = int)
parser.add_argument(
    "--min_tsse", 
    default = 1,
    type = int)
parser.add_argument(
    "--set", 
    action='store_true',
    help='whether read from sets of adata')

genome = {'mm10': snap.genome.mm10, 'hg38': snap.genome.hg38}
anno_gff = {
    "mm10": os.environ.get("MM10_GFF3", "/path/to/gencode.mm10.annotation.gff3"),
    "hg38": os.environ.get("HG38_GFF3", "/path/to/gencode.hg38.annotation.gff3"),
}


def label_transfer(atac_st, reference, query, output_dir, sampleid):
    query.obs['cell_type'] = pd.NA
    data = ad.concat(
        [reference, query],
        join='inner',
        label='batch',
        keys=["reference", "query"],
        index_unique='_'
    )

    sc.pp.filter_genes(data, min_cells=5)
    sc.pp.highly_variable_genes(
        data,
        n_top_genes = 2000,
        flavor="seurat_v3",
        batch_key="batch",
        subset=True
    )
    lvae_path = os.path.join(output_dir, sampleid + "_saved_model")
    if not os.path.exists(lvae_path):
        print(f'Modeling trianing...', file=sys.stderr)
        scvi.model.SCVI.setup_anndata(data, batch_key="batch")
        vae = scvi.model.SCVI(
            data,
            n_layers=2,
            n_latent=15, # Dimensionality of the latent space
            gene_likelihood="nb",
            dispersion="gene-batch",
        )

        vae.train(max_epochs=1000, early_stopping=True)

        # vae.history['elbo_validation'].plot(ax=ax)

        temp = reference.obs
        temp.index = temp.index + "_reference"
        data.obs["celltype_scanvi"] = 'Unknown'
        ref_idx = data.obs['batch'] == "reference"
        data.obs["celltype_scanvi"][ref_idx] = temp['cell_type'][ref_idx]

        lvae = scvi.model.SCANVI.from_scvi_model(
            vae,
            adata=data,
            labels_key="celltype_scanvi",
            unlabeled_category="Unknown",
        )

        lvae.train(max_epochs=1000, n_samples_per_label=100)

        lvae.save(lvae_path, save_anndata=True)

        # lvae.history['elbo_train'][1:].plot()
    else:
        print(f'Modeling training finished, loading {lvae_path}...', file=sys.stderr)
        lvae = scvi.model.SCVI.load(lvae_path)

    data.obs["C_scANVI"] = lvae.predict(data)
    data.obsm["X_scANVI"] = lvae.get_latent_representation(data)

    sc.pp.neighbors(data, use_rep="X_scANVI")
    sc.tl.umap(data)
    
    new_name = [i + "_query" for i in atac_st.obs_names]
    atac_st.obs['cell_type'] = data.obs.loc[new_name]['C_scANVI'].to_numpy() 

    # atac_st.obs['cell_type'] = data.obs.loc[atac_st.obs_names + '_query']['C_scANVI'].to_numpy() 
    
    return atac_st, data


def main():
    args = parser.parse_args()

    print(f'input_h5ad: {args.input_h5ad}', file=sys.stderr)
    print(f'reference_h5ad: {args.reference_h5ad}', file=sys.stderr)
    print(f'build: {args.species}', file=sys.stderr)  
    print(f'gene_anno: {anno_gff[args.species]}', file=sys.stderr)
    if args.input_h5ad[-4::] != "h5ad":
        flag = 0
        print('Importing fragment file...', file=sys.stderr)

        atac_st = snap.pp.import_data(
                args.input_h5ad,
                chrom_sizes=genome[args.species],
                sorted_by_barcode=False,
                min_num_fragments=10,
                n_jobs=-1
            )
        snap.metrics.tsse(atac_st, gene_anno = anno_gff[args.species])

        out_tss = args.input_h5ad.replace("_fragments.tsv.gz", ".tss_enrichment.pdf")
        snap.pl.tsse(atac_st,
            interactive=False, 
            out_file=out_tss, 
            min_fragment=10)
        out_metadata = args.input_h5ad.replace("_fragments.tsv.gz", ".metadata.txt")

        atac_st.obs.to_csv(out_metadata, sep="\t",header=True)

        ## clustering
        snap.pp.filter_cells(atac_st, 
            min_counts=args.min_fragments_n,
            min_tsse=args.min_tsse,
            max_counts=args.max_fragments_n)
        snap.pp.add_tile_matrix(atac_st)
        snap.pp.select_features(atac_st, n_features=250000, max_iter = 2)
        
        snap.tl.spectral(atac_st)
        snap.tl.umap(atac_st)
        snap.pp.knn(atac_st)
        snap.tl.leiden(atac_st)
        outpdf = args.input_h5ad.replace("_fragments.tsv.gz", ".cluster.pdf")
        snap.pl.umap(atac_st, color='leiden', interactive=False, out_file=outpdf, height=500)
    else:
        flag = 1
        if args.set:
            atac_st = snap.read_dataset(args.input_h5ad)
        else:
            atac_st = snap.read(args.input_h5ad, backed=None)
    reference = sc.read_h5ad(args.reference_h5ad)

    if args.min_fragments_n != 1000 and flag == 1:
        snap.pp.filter_cells(atac_st,
            min_counts=args.min_fragments_n,
            min_tsse=args.min_tsse,
            max_counts=args.max_fragments_n)

        snap.pp.add_tile_matrix(atac_st)
        snap.pp.select_features(atac_st, n_features=250000, max_iter = 2)
        snap.tl.spectral(atac_st)
        snap.tl.umap(atac_st)
        snap.pp.knn(atac_st)
        snap.tl.leiden(atac_st)

    query = snap.pp.make_gene_matrix(atac_st, anno_gff[args.species])
    atac_st_annotated, data = label_transfer(atac_st, 
        reference,
        query,
        args.output_dir,
        args.sampleid )
    atac_st.close()
    
    if args.input_h5ad[-4::] != "h5ad":
        out_file_h5ad = args.input_h5ad.replace("_fragments.tsv.gz", ".annotated.h5ad")
        out_coemedding_h5ad = args.input_h5ad.replace("_fragments.tsv.gz", ".co_embedding.h5ad")

    else:
        out_file_h5ad = args.input_h5ad.replace(".h5ad", ".annotated.h5ad")
        out_coemedding_h5ad = args.input_h5ad.replace(".h5ad", ".co_embedding.h5ad")

    atac_st_annotated.write(out_file_h5ad, compression="gzip")
    data.write(out_coemedding_h5ad, compression="gzip")

    output_prefix = "_" + args.sampleid
    sc.pl.umap(data,
        color=['C_scANVI', "batch"],
        wspace=0.45,
        save = output_prefix + ".coembedding.pdf")

    sc.pl.umap(atac_st_annotated, 
        color=['cell_type', "leiden"],
        wspace=0.45, 
        save = output_prefix + ".celltype_atac.pdf")

    gene_matrix = snap.pp.make_gene_matrix(atac_st_annotated, anno_gff[args.species])


    sc.pp.filter_genes(gene_matrix, min_cells= 5)
    sc.pp.normalize_total(gene_matrix)
    sc.pp.log1p(gene_matrix)
    sc.external.pp.magic(gene_matrix, solver="approximate")

    gene_matrix.obsm["X_umap"] = atac_st_annotated.obsm["X_umap"]
    marker_genes = [
        'Slc17a7', # Glutamatergic neuron
        'Gad2', # GABAergic neuron         
        'Arpp21', 'Nrg3', # Excitatory Neuron
        'Snhg11', 'Pvalb', 'Lamp5', 'Vip', 'Grip1', # Inhibitory Neuron
        'Slc1a3', 'Slc1a2','Aldh1l1', # Astrocyte
        'Rgs9', 'Rarb', # Medium_spiny_neurons
        'Tgfbr1', 'Srgap2', 'P2ry12', # Microglial_cells
        'Lhfpl3', 'Dscam', # Oligodendrocyte_precursor_cells
        'Mbp', 'Plp1', 'Mog' # Oligodendrocyte
                ]
    sc.pl.umap(gene_matrix, 
        use_raw=False, 
        color= ["cell_type"] + marker_genes,
        save = output_prefix + ".marker_gene.pdf" )


if __name__ == '__main__':
    main()

