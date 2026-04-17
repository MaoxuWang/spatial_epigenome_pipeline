suppressPackageStartupMessages({
    library(Seurat)
    library(trqwe)
    library(BuenColors)
    library(tidyverse)
    library(Signac)
    library(rhdf5)
    library(SingleCellExperiment)
    library(chromVAR)
    library(motifmatchr)
    library(chromVARmotifs)
    library(BiocParallel)
    library(GenomeInfoDb)
})

args = commandArgs(T)

peak_h5_path = args[1]
fragment_path = args[2]
species = args[3]

read_h5 <- function(
    h5_file){

    # read sparse matrix (dgCMatrix)
    h5_file <- H5Fopen(h5_file)
    counts <- as(h5_file$mat$block0_values, "dgCMatrix")  
    rownames(counts) <- h5_file$mat$axis0  # "chr1:100-200"
    colnames(counts) <- h5_file$mat$axis1
    H5Fclose(h5_file)

    return(counts)
} 


createSignac <- function(
    counts,
    fragment_path,
    species
){
    if (species == "mm10"){
        chrom_info <- read.table("~/database/chromosomeSize/mm10.chrom.sizes", header = FALSE)
        # create seqinfo
        seqinfo <- Seqinfo(seqnames = chrom_info$V1, seqlengths = chrom_info$V2, genome = "mm10")

        library(EnsDb.Mmusculus.v79)
        annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
        seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
        genome(annotations) <- "mm10"
    }else{
        chrom_info <- read.table("~/database/chromosomeSize/hg38.autosome.chrom.sizes", header = FALSE)
        # create seqinfo
        seqinfo <- Seqinfo(seqnames = chrom_info$V1, seqlengths = chrom_info$V2, genome = "hg38")

        library(EnsDb.Hsapiens.v86)
        annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
        seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
        genome(annotations) <- "hg38"
    }

    signac_assay <- CreateChromatinAssay(
        counts = counts,
        sep = c(":", "-"),  # separate
        genome = seqinfo,
        fragments = fragment_path,
        min.cells = 50
    )

    signac_obj <- CreateSeuratObject(
        counts = signac_assay,
        assay = 'peaks',
        project = 'ATAC'
    )

    # add the gene information to the object
    Annotation(signac_obj) <- annotations
    signac_obj %>%
        FindTopFeatures(min.cutoff = 50) %>%
        RunTFIDF() -> signac_obj

    signac_obj %>%
        RunSVD() %>%
        RunUMAP(reduction = 'lsi', dims = 1:15, verbose = FALSE) %>%
        FindNeighbors(reduction = 'lsi', dims = 1:15, verbose = FALSE) %>%
        FindClusters(algorithm = 3, resolution = c(0.2, 0.4), verbose = FALSE) -> signac_obj
    
    extra_r_dir <- Sys.getenv("DSCHIC_R_FUNCTION_DIR", unset = "")
    if (nzchar(extra_r_dir)) {
        source(file.path(extra_r_dir, "scRNA.func.R"))
        source(file.path(extra_r_dir, "plot_function.R"))
    }

    p <- plot_seurat(signac_obj,
        signac_obj$seurat_clusters,
        "seurat_clusters",
        "Peak clustering",
        label=FALSE,
        raster = FALSE,
        point_size = .6)
    ggsave(
        p, 
        filename = gsub(".fragments.tsv.gz", ".peak_clustering.pdf", fragment_path),
        units = "mm", width = 120, height = 120, dpi = 600
    )

    return(signac_obj)
}


process_chromVAR <- function(
    counts,
    peak_h5_path,
    species
){
    if(species == "mm10"){
        library(BSgenome.Mmusculus.UCSC.mm10)
        genome = BSgenome.Mmusculus.UCSC.mm10
    }else{
        library(BSgenome.Hsapiens.UCSC.hg38)
        genome = BSgenome.Hsapiens.UCSC.hg38
    }
    peaks = StringToGRanges(rownames(counts), sep = c(":", "-"))
    fragment_counts <- SummarizedExperiment(
        assays = list(counts = counts),
        rowRanges = peaks
    )   
    BiocParallel::register(BiocParallel::MulticoreParam(8, progressbar = TRUE))

    fragment_counts <- addGCBias(
        fragment_counts,
        genome=genome)
    mcsaveRDS(fragment_counts,
        file = gsub(".peak_matrix.h5", ".fragment_counts.rds", peak_h5_path))

    # Get background peaks for chromVAR
    bg <- getBackgroundPeaks(
        object = fragment_counts,
        niterations = 200)

    # K-mer scoring (used for clustering single cells)
    kmer_ix <- matchKmers(
        k = 6,
        subject = fragment_counts,
        genome=genome
    )
    # Motif scoring (used to annotate single cells)
    motif_ix <- matchMotifs(
        pwms = mouse_pwms_v2,
        subject = fragment_counts,
        genome=genome
    )
    print("Computing deviations for motifs ..\n\n")
    dev_motif <- computeDeviations(object = fragment_counts,annotations = motif_ix,background_peaks=bg)

    print("Computing deviations for k-mers ...")
    dev_kmer <- computeDeviations(object = fragment_counts,annotations = kmer_ix,background_peaks=bg)

    kmer_6_activity <- deviationScores(dev_kmer)
    mcsaveRDS(kmer_6_activity, 
        file = gsub(".peak_matrix.h5", ".kmer_6_activity.rds", peak_h5_path))

    motif_activity <- deviationScores(dev_motif)
    mcsaveRDS(motif_activity, 
        file = gsub(".peak_matrix.h5", ".motif_activity.rds", peak_h5_path))
}


counts = read_h5(peak_h5_path)
mcsaveRDS(counts, 
    file = gsub(".peak_matrix.h5", ".peak_matrix.rds", peak_h5_path))

signac_obj <- createSignac(
    counts,
    fragment_path,
    species
)
mcsaveRDS(signac_obj, 
    file = gsub(".peak_matrix.h5", ".signac.rds", peak_h5_path))

process_chromVAR(
    counts,
    peak_h5_path,
    species
)