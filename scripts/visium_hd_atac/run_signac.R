library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79) # Mouse brain genome annotation
library(argparse)
library(tidyverse)
library(dplyr)
library(GenomicRanges)
library(GenomicAlignments)
library(Rsamtools)
library(trqwe)
library(scales)
library(pals)
library(future) # For parallel processing
library(doParallel) # For parallel processing
library(foreach) # For parallel processing
library(BiocParallel)
library(data.table)

cat("--- R script starting: Signac spatial epigenomics pipeline (SE IVT deduplication version) ---\n")

parser <- ArgumentParser()
parser$add_argument("--input_file", type="character", required=TRUE, 
                    help="[Required] Raw SE BAM file with barcodes (e.g., possorted_genome_bam.bam)")
parser$add_argument("--spatial_dir", type="character", required=TRUE, 
                    help="[Required] Path to Visium-HD 'spatial' directory (contains tissue_positions.parquet and image files)")
parser$add_argument("--metadata_file", type="character", 
                    help="[Optional] Tab-separated file with additional metadata (row names as barcodes)")
parser$add_argument("--outdir", type="character", required=TRUE, 
                    help="[Required] Output directory for saving .rds and .pdf files")
parser$add_argument("--broad", type="logical", default=FALSE, 
                    help="[Required] Output directory for saving .rds and .pdf files")
parser$add_argument("--sample_name", type="character", required=TRUE, 
                    help="[Required] Sample name (e.g., 'my_sample')")
parser$add_argument("--genome", type="character", default="mm10", 
                    help="Genome version (e.g., 'mm10')")
parser$add_argument("--extend_bp", type="integer", default=200, 
                    help="Base pairs to extend SE reads; ~200bp is biologically appropriate for H3K4me3")
parser$add_argument("--tss_upstream", type="integer", default=2000, 
                    help="Distance upstream of TSS for Gene Activity calculation")
parser$add_argument("--mapq_filter", type="integer", default=30, 
                    help="min mapping quality required")
parser$add_argument("--resolution", type="double", default=0.4, 
                    help="Resolution for downstream clustering")
parser$add_argument("--tag", type="character", default="CB", 
                    help="BAM tag for cell barcode (e.g., CB, CR, XC)")
parser$add_argument("--n_cores", type="integer", default=1, 
                    help="Number of cores to use for parallel processing")
parser$add_argument("--assay", type="character", default="H3K4me3", 
                    help="CUTTAG or ATAC peak matrix assay")
parser$add_argument("--decontam_path", type="character", 
                    default=Sys.getenv("REMOVE_SPOT_SWAP", "tools/remove_spot_swap.sh"), 
                    help="script to remove contamination in spot")

Args <- parser$parse_args()

effective_size <- c("hg38" = 2913022398, "mm10" = 2652783500)


readImage <- function(
    img_dir, assay){
    position_file = paste0(img_dir, "/tissue_positions.parquet")
    

    if (file.exists(position_file)) {
        df <- read_parquet(position_file)

        df %>%
            mutate(barcode = gsub("-1", "", barcode)) %>%
            write.table(file = paste0(img_dir, "/tissue_positions_list.csv"),
                sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
        # 执行重命名
        file.rename(from = position_file, to = paste0(img_dir, "/tissue.cr.positions.parquet"))
    } else {
        print("parquet file is already processed!")
    }   
        
    image <- Read10X_Image(
        image.dir = img_dir, 
        assay = assay,
        image.type = "VisiumV2",
        filter.matrix = TRUE)

    return(image)
}


cat(paste("Setting up parallel processing with", Args$n_cores, "cores...\n"))
if (Args$n_cores > 1) {
  plan("multicore", workers = Args$n_cores)
  options(future.globals.maxSize = 50000 * 1024^2) #50GB memory limit
  register(MulticoreParam(workers = Args$n_cores))
  registerDoParallel(cores = Args$n_cores)
}

dir.create(Args$outdir, recursive = TRUE, showWarnings = FALSE)
output_rds <- file.path(Args$outdir, paste0(Args$sample_name, ".signac.rds"))
output_raw_rds <- file.path(Args$outdir, paste0(Args$sample_name, ".signac_raw.rds"))

output_pdf <- file.path(Args$outdir, paste0(Args$sample_name, ".n_fragment_qc.pdf"))


cat(paste("Input BAM:", Args$input_file, "\n"))
cat(paste("Spatial directory:", Args$spatial_dir, "\n"))
cat(paste("Output RDS:", output_rds, "\n"))
cat(paste("Using barcode tag:", Args$tag, "\n"))

## prepare annotation
if (Args$genome == "mm10"){
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

if (grepl("\\.bam$", Args$input_file, ignore.case = TRUE)) {
  output_ext_frag_path <- file.path(Args$outdir, paste0(Args$sample_name, ".extended_fragments.tsv.gz"))
  output_ext_6col_frag_path <- file.path(Args$outdir, paste0(Args$sample_name, ".extended_fragments_6_col.tsv"))

  output_tsv_path <- file.path(Args$outdir, paste0(Args$sample_name, ".extended_fragments.tsv"))

  if (file.exists(output_ext_frag_path)) {
    
    cat(paste("--- Found existing extended file: ", output_ext_frag_path, "---\n"))
    
    fragments_df <- data.table::fread(output_ext_frag_path, nThread = Args$n_cores)
    colnames(fragments_df) <- c("chr", "start", "end", "barcode", "count")
    extended_fragments <- makeGRangesFromDataFrame(
      fragments_df,
      keep.extra.columns = TRUE,
      starts.in.df.are.0based = TRUE # <-- Critical step
    )
    rm(fragments_df)


  } else {
      
    cat(paste0("Filtering mapping quality ", Args$mapq_filter, " and dedup bam...\n"))

    bam_file_path <- Args$input_file
    bf <- BamFile(bam_file_path)
    seq_info <- seqinfo(bf)

    chromosomes_to_process <- seqnames(seq_info)[grepl("^(chr[0-9XYM]+$)", seqnames(seq_info))]
    cat(paste("将处理以下染色体:", paste(chromosomes_to_process, collapse=", "), "\n"))

    base_param <- ScanBamParam(
        what = "seq", 
        tag = Args$tag,
        mapqFilter = Args$mapq_filter,
        flag = scanBamFlag(
            isSecondaryAlignment = FALSE, 
            isUnmappedQuery = FALSE,
            isNotPassingQualityControls = FALSE
        )
    ) 

    # 使用 foreach 按染色体并行处理，并用 rbind 合并
    gr_counts_df <- foreach(chr = chromosomes_to_process, .combine = 'rbind', .packages = c("GenomicRanges", "GenomicAlignments", "Rsamtools", "dplyr", "tibble")) %dopar% {
        cat(paste("Processing ", chr, "...\n"))
        chr_param <- base_param
        which_gr <- GRanges(seqnames = chr, ranges = IRanges(1, seqlengths(seq_info)[chr]))
        bamWhich(chr_param) <- which_gr
        
        aln_chr <- readGAlignments(bam_file_path, param = chr_param)
        gr_chr <- GRanges(aln_chr)
        
        if (is.null(mcols(gr_chr)[[Args$tag]])) {
            cat(paste("  - WARNING: ", chr, " not found '", Args$tag, "' tags\n"))
            return(NULL)
        }

        # 转换为 Tibble
        gr_df_chr <- tibble(
          seqnames = as.character(seqnames(gr_chr)),
          start = start(gr_chr),
          end = end(gr_chr), 
          strand = as.character(strand(gr_chr)),
          barcode = mcols(gr_chr)[[Args$tag]]
        )
        
        gr_df_chr_with_5prime <- gr_df_chr %>%
          mutate(five_prime_coord = ifelse(strand == "+", start, end))

        gr_counts_df_chr <- gr_df_chr_with_5prime %>%
          group_by(barcode, seqnames, five_prime_coord, strand) %>%
          summarise(count = n(), .groups = 'drop') %>%
          mutate(start = five_prime_coord, end = five_prime_coord) %>% # <-- [CHANGED]
          dplyr::select(seqnames, start, end, barcode, count, strand)
        
        return(gr_counts_df_chr)
    }

    cat("\nBAM processed successfully!\n")
    cat(paste("  - Number of unique fragments:", nrow(gr_counts_df), "\n"))

    cat("Step 3: extend single-end fragments", Args$extend_bp, "bp...\n")

    gr_dedup <- makeGRangesFromDataFrame(
      gr_counts_df,
      seqnames.field = "seqnames",
      start.field = "start",
      end.field = "end",
      strand.field = "strand",
      keep.extra.columns = TRUE, # keep.extra.columns 会保留 'barcode' 和 'count' 列
      starts.in.df.are.0based = FALSE
    )
    
    extended_fragments <- Signac::Extend(
      gr_dedup,
      upstream = 0,
      downstream = Args$extend_bp
    )
    # mcols(extended_fragments)$barcode 已经存在
    rm(gr_counts_df, gr_dedup)
    

    cat("...Converting extended fragments to 0-based BED6 format...\n")
    fragments_df <- data.table(
        chr = as.character(seqnames(extended_fragments)),
        start = sprintf("%.0f", start(extended_fragments) - 1), # *** 0-BASED START ***
        end = sprintf("%.0f", end(extended_fragments)),
        barcode = extended_fragments$barcode,
        count = extended_fragments$count,
        strand = as.character(strand(extended_fragments))
    )

    chrs_in_data <- unique(fragments_df$chr)
    
    # 按照 'chromosomes_to_process' 的顺序来过滤和排序
    ordered_chr_levels <- intersect(chromosomes_to_process, chrs_in_data)

    fragments_df[, chr := factor(chr, levels = ordered_chr_levels)]

    cat("...Sorting fragments for tabix...\n")
    if (!is.numeric(fragments_df$start)) {
      fragments_df[, start := as.numeric(as.character(start))]
      fragments_df[, end := as.numeric(as.character(end))]

    }

    setkey(fragments_df, chr, start) # Sort by chr and start
    options(scipen = 999)
    cat(paste("...Saving fragments to:", output_ext_frag_path, "...\n"))
    
    data.table::fwrite(
      fragments_df, 
      output_ext_6col_frag_path, 
      sep="\t", 
      compress="none", 
      nThread = Args$n_cores, 
      col.names = FALSE
    )

    cat("...Save complete.\n")
    fragments_df_signac <- fragments_df %>% dplyr::select(-strand)
    data.table::fwrite(
      fragments_df_signac, 
      output_tsv_path, 
      sep="\t", 
      compress="none", 
      nThread = Args$n_cores, 
      col.names = FALSE
    )

    Rsamtools::bgzip(file = output_tsv_path, dest = output_ext_frag_path)

    # *** NEW: Create Tabix index ***
    cat("...Indexing fragments file with Rsamtools::indexTabix...\n")
    Rsamtools::indexTabix(file = output_ext_frag_path, format = "bed")
    cat("...Tabix index created.\n")
    
    # We no longer need the GRanges object
    rm(fragments_df, extended_fragments) 
    file.remove(output_tsv_path) # (可选) 删除临时的未压缩文件
    cat("...save successfully\n")
  }

  cat("--- [Integration Step] Running AWK script for cross-contamination removal... ---\n")

  # 1. 定义清理后的文件路径
  decontam_frag_path <- sub("\\.tsv\\.gz$", ".decontaminated.fragments.tsv.gz", output_ext_frag_path)
  # (如果输入是 .fragments.tsv.gz, 这会创建一个 .decontaminated.fragments.tsv.gz)
  decontam_script_path <- Args$decontam_path
  # 2. 定义 decontamination 脚本的路径 (假设它在 'scripts' 目录中)
  #    *** 重要: 请根据需要修改此路径 ***
  if (!file.exists(decontam_script_path)) {
    stop(paste("Decontamination script not found at:", decontam_script_path))
  }

  # 3. 构建并执行 system() 命令
  command <- paste("bash", decontam_script_path, output_ext_6col_frag_path, decontam_frag_path)

  cat(paste("Executing:", command, "\n"))
  return_code <- system(command) # <-- 执行 bash 脚本

  if (return_code != 0) {
      stop(paste("Decontamination script failed with exit code:", return_code))
  }

  cat(paste("--- Contamination removal complete. Cleaned file:", decontam_frag_path, "---\n"))

  # 4. *** 关键: 将主 fragments 路径更新为清理后的文件 ***
  output_ext_frag_path <- decontam_frag_path

} else if (grepl("\\.fragments\\.tsv\\.gz$", Args$input_file, ignore.case = TRUE)) {

    cat("--- Input is a .fragments.tsv.gz file. Skipping Steps 1-3. ---\n")
    
    # 将 output_ext_frag_path 直接设置为输入的 fragments 文件
    output_ext_frag_path <- Args$input_file
    
    # *** 关键检查 ***: 确保该文件已被 tabix 索引
    tbi_file <- paste0(output_ext_frag_path, ".tbi")
    if (!file.exists(tbi_file)) {
        cat(paste("--- ERROR: Index file not found:", tbi_file, "---\n"))
        cat("--- Please index your fragments file before proceeding. ---\n")
        cat("--- You can use: Rsamtools::indexTabix(file = \"", output_ext_frag_path, "\", format = \"bed\") ---\n")
        stop("Input fragments file must be tabix-indexed (.tbi file missing)!")
    }
    
    cat(paste("--- Using existing, indexed fragments file:", output_ext_frag_path, "---\n"))

} else {
    stop(paste("Unknown input file type:", Args$input_file, 
               ". Please provide a .bam or .fragments.tsv.gz file."))
}

cat("Step 4: Creating ChromatinAssay...\n")

print(paste("Fragment file used in downstream analysis: ", output_ext_frag_path))
## call peak first
fragments_obj <- CreateFragmentObject(
  path = output_ext_frag_path 
)
cat(paste("Call peak (Broad set to be: ", Args$broad, ")\n"))

peaks <- CallPeaks(
  object = fragments_obj,
  genome = Args$genome,
  outdir = Args$outdir,
  effective.genome.size = effective_size[Args$genome],
  extsize = Args$extend_bp,
  shift = 0,
  broad = Args$broad,
  additional.args = "--nomodel --keep-dup all --qvalue 0.01",
  macs2.path = Sys.getenv("MACS2", "macs2")  # 确保MACS2已安装并指定正确路径
)

## create peak matrix

counts_matrix <- FeatureMatrix(
  fragments = fragments_obj,
  features = peaks,
  process_n = 10000 
)

chrom_assay <- CreateChromatinAssay(
  counts = counts_matrix,
  genome = seqinfo,
  fragments = output_ext_frag_path,
  min.cells = 100,
  min.features = 0)

cat("Step 5: Creating Seurat object...\n")
sp_obj <- CreateSeuratObject(
  counts = chrom_assay, # Create an empty object
  project = Args$sample_name,
  assay = Args$assay
)

cat("Step 6: Loading spatial data...\n")
image <- readImage(Args$spatial_dir, Args$assay)

DefaultAssay(object = image) <- Args$assay

cell_used <- intersect(Cells(image), Cells(sp_obj))
image <- image[cell_used]
sp_obj <- subset(sp_obj, cells = cell_used)

sp_obj[['slice_1']] <- image
Annotation(sp_obj) <- annotations

cat("Step 7: Calculating TSS enrichment and FRiP...\n")
sp_obj <- PercentageFeatureSet(sp_obj, pattern = "^chrM-", col.name = "percent.mt", assay = Args$assay)
sp_obj <- TSSEnrichment(sp_obj, assay = Args$assay)

if (!is.null(Args$metadata_file)) {
  cat("Step 8: Adding external metadata...\n")
  metadata <- read.table(Args$metadata_file, 
                         sep = "\t", 
                         row.names = 1, 
                         header = TRUE)
  sp_obj <- AddMetaData(sp_obj, metadata = metadata)
}
sp_obj <- FRiP(
  object = sp_obj,
  assay = Args$assay,
  total.fragments = 'n_fragment'
)
cat("Step 9: Calculating Gene Activity Score...\n")
# GeneActivity will use our provided, already deduplicated and extended "fragments"
gene_activities <- GeneActivity(
  object = sp_obj,
  assay = Args$assay,
  extend.upstream = Args$tss_upstream 
)

# Add Gene Activity Score as a new "RNA" Assay
sp_obj[['RNA']] <- CreateAssayObject(counts = gene_activities)

cat("Step 10: Filtering and normalizing Gene Activity (RNA Assay)...\n")
valid_genes <- grep("^Gm|^Mir", rownames(sp_obj[['RNA']]), invert = TRUE, value = TRUE)
sp_obj[['RNA']] <- subset(sp_obj[['RNA']], features = valid_genes)
mcsaveRDS(sp_obj, output_raw_rds) # Save final Seurat object

DefaultAssay(sp_obj) <- "RNA"
cat("Step 11: Running PCA and clustering (based on RNA Assay)...\n")
plan("sequential")

sp_obj %>%
    NormalizeData(verbose = FALSE, assay = 'RNA', normalization.method = 'LogNormalize') %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(verbose = FALSE, dims = 1:10) %>%
    FindClusters(verbose = FALSE, resolution = Args$resolution) %>%
    RunUMAP(verbose = FALSE, dims = 1:10) -> sp_obj

cat("Step 12: Saving RDS object and QC plots...\n")
mcsaveRDS(sp_obj, output_rds) # Save final Seurat object

data_range <- quantile(sp_obj$nCount_H3K4me3, probs = c(0.02, 0.98), na.rm = TRUE)

p <- SpatialFeaturePlot(sp_obj,
    features = "nCount_H3K4me3",
    image.alpha = 0, 
    stroke = 0, 
    max.cutoff = "q99",
    min.cutoff = "q5",
    pt.size.factor = 2.5)

p <- p & scale_fill_gradientn(
  name = "Fragment Count",
  colors = inferno(100),

  # We still need to specify where to place ticks
  breaks = scales::pretty_breaks(n = 3),
    
  # labels parameter is now a function that automatically handles values in breaks
  labels = label_number(scale = 1e-3, suffix = "K"),
  oob = scales::squish,
  limits = data_range
  
) & theme(
  legend.title = element_text(face = "bold", size = 12, hjust = 1.5),
  legend.text = element_text(size = 10)
)
p
ggsave(p, filename = output_pdf, 
       width = 120, height = 150, units = "mm",
       create.dir = TRUE, dpi = 300, limitsize = FALSE)

cat("--- R script execution completed ---\n")