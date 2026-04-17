suppressPackageStartupMessages({library(Seurat)
library(tidyverse)
library(trqwe)
library(argparse)})

options(future.globals.maxSize = 8000 * 1024^2)

parser <- ArgumentParser()
parser$add_argument("--count_matrix_path", type="character", help="count matrix file")
parser$add_argument("--metadata_file", type="character", help="metadata file")
parser$add_argument("--spatial_barcode_file", type="character", help="spatial barcode file")
parser$add_argument("--out_seurat", type="character", help="output seurat object file")
parser$add_argument("--image_dir", type="character", help="image directory")
parser$add_argument("--in_tissue", type="logical", help="whether in tissue")
Args <- parser$parse_args()

in_tissue = as.logical(Args$in_tissue)

counts = read_delim(Args$count_matrix_path, delim = "\t")

counts <- data.frame(counts, check.names = FALSE)
colnames(counts)[1] = "index"

valid_genes <- grep("^Gm|^Mir", 
                    colnames(counts), 
                    invert = TRUE, value = TRUE)

counts <- counts[, c("index", valid_genes)]
print(paste0("Gene used: ", length(valid_genes)))

sampleid = gsub(".metadata.txt", "", gsub(".*\\/", "", Args$metadata_file))
spatial_barcode <- read_delim(Args$spatial_barcode_file, delim = "\t")
colnames(spatial_barcode) <- c("cellbarcode", "index")

counts %>%
    inner_join(spatial_barcode, by = "index") %>%
    column_to_rownames("cellbarcode") %>%
    dplyr::select(-index) -> counts


metadata  <- read.table(Args$metadata_file, 
    sep = "\t", 
    row.names = 1, 
    header = TRUE)
metadata$cellbarcode = rownames(metadata)

createSpatialObj <- function(
    counts,
    image_dir,
    in_tissue,
    metadata,
    out_seurat){

    print("Reading image...")
    print(paste0("In tissue: ", in_tissue))

    image <- Read10X_Image(
        image.dir = image_dir,
        image.type = "VisiumV2",
        filter.matrix = in_tissue)
    image@assay <- "RNA"

    counts <- counts[Cells(image), ]
    object_spatial <- CreateSeuratObject(t(counts), 
        min.cells = 50,
        project = sampleid)

    object_spatial <- AddMetaData(object_spatial, metadata = metadata)

    object_spatial[['slice_1']] <- image
    p <- SpatialFeaturePlot(object_spatial,
        features = "n_fragment",
        image.alpha = 0, 
        stroke = 0.2, 
        pt.size.factor = 3) + 
        theme(
            legend.title = element_blank(), 
            legend.text = element_text(size=12)) +  # 设置legend标签之间的大小
        guides(color = guide_legend(override.aes = list(size=5)))

    p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)

    pc.num = 1:20

    # object_spatial %>%
    #     SCTransform(assay = "RNA", return.only.var.genes = FALSE, 
    #         vars.to.regress = "nCount_RNA",
    #         verbose = FALSE ) %>%
    #     RunPCA(assay = "SCT", verbose = FALSE, pcs.compute=50) %>%
    #     FindNeighbors(dims = 1:25, verbose = FALSE) %>%
    #     FindClusters(verbose = FALSE, resolution = c(0.2, 0.4,0.6)) %>%
    #     RunUMAP(dims = 1:25, verbose = FALSE) -> object_spatial
    if (in_tissue == TRUE){
        out_n_frag_pdf = gsub(".rds", ".n_frag.pdf", out_seurat)

        object_spatial %>%
            NormalizeData(verbose = FALSE) %>%
            FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
            # ScaleData(verbose = FALSE, vars.to.regress = "nCount_RNA") %>%
            ScaleData(verbose = FALSE) %>%
            RunPCA(verbose = FALSE) %>%
            RunUMAP(reduction = "pca", dims = pc.num, verbose = F) %>%
            FindNeighbors(reduction = "pca", dims = pc.num, verbose = F) %>% 
            FindClusters(resolution=c(0.2, 0.4, 0.6)) -> object_spatial
    }else{
        out_n_frag_pdf = gsub(".rds", ".all_n_frag.pdf", out_seurat)

    }
     ggsave(p, filename = out_n_frag_pdf, units = "mm",
            width = 120, height = 120, dpi = 600)
    return(object_spatial)
}

object_spatial.in_tissue <- createSpatialObj(
    counts,
    Args$image_dir,
    in_tissue,
    metadata,
    Args$out_seurat
    )
mcsaveRDS(object_spatial.in_tissue, Args$out_seurat)