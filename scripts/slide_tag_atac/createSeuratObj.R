library(SeuratDisk)
library(Seurat)
library(tidyverse)
library(trqwe)
library(arrow)
library(argparse)


parser <- ArgumentParser()
parser$add_argument("--input_h5ad", required = TRUE,
    help = "gene_activity.h5ad")
parser$add_argument("--metadata_file", required = TRUE,
    help = "metadata.txt")
parser$add_argument("--position_file", required = TRUE,
    help = "cell_locations.tsv.gz")
parser$add_argument("--img", required = TRUE,
    help = "_aligned_DAPI.png")
parser$add_argument("--size_x", default = 0,
    help = "Required. Image size x")
parser$add_argument("--size_y", default = 0,
    help = "Required. Image size y")
parser$add_argument("--outdir", default = "./",
    help = "The output directory")

Args <- parser$parse_args()


sampleid = gsub(".gene_activity.h5ad", "", Args$input_h5ad) # out/sample_name
sampleName = gsub(".*\\/", "", sampleid)
sampleid = paste0(Args$outdir,"/", sampleName)

## Step 1. convert h5ad to seurat

Convert(Args$input_h5ad, 
        meta.data = FALSE,
        misc = FALSE,
        dest = "h5seurat", overwrite = TRUE)

input_h5seurat <- gsub("h5ad", "h5seurat", Args$input_h5ad)

object_spatial <- LoadH5Seurat(input_h5seurat, meta.data = F)
valid_genes <- grep("^Gm|^Mir", 
                    rownames(object_spatial), 
                    invert = TRUE, value = TRUE)

object_spatial <- subset(object_spatial, features = valid_genes)

metadata <- read.table(Args$metadata_file, 
    sep = "\t", 
    row.names = 1, 
    header = TRUE)
metadata$cellbarcode = rownames(metadata)

object_spatial@meta.data %>%
    mutate(cellbarcode = rownames(object_spatial@meta.data)) %>%
    inner_join(metadata) -> object_spatial@meta.data 
rownames(object_spatial@meta.data) <- object_spatial@meta.data$cellbarcode

## Step 2. make spatial object
df_center = read.csv(Args$position_file, sep = "\t")
df_center = df_center %>% filter(Cell_Barcode %in% Cells(object_spatial))
barcodes = df_center$Cell_Barcode
spatial_x = df_center$X
spatial_y = df_center$Y

object_spatial <- subset(object_spatial, cells = barcodes)

### create spatial reduction
spatial_embeddings <- matrix(data = c(spatial_x, spatial_y), nrow = length(spatial_x), ncol = 2, byrow = FALSE)
rows = barcodes
cols = c("spatial_1", "spatial_2")
attr(spatial_embeddings, "dimnames") <- list(rows, cols)

dim_reduc_obj <- CreateDimReducObject(embeddings = spatial_embeddings, assay = "RNA", key = "spatial")
object_spatial@reductions$spatial = dim_reduc_obj

### sort
spatial_matrix <- object_spatial@reductions$spatial@cell.embeddings
spatial_matrix_sorted <- spatial_matrix[match(row.names(object_spatial@meta.data), row.names(spatial_matrix)),]
object_spatial@reductions$spatial<- CreateDimReducObject(embeddings = spatial_matrix_sorted, key='spatial_', assay='RNA')

### add image layer
samplename = 'slice_1'
img_64 = base64enc::dataURI(file = Args$img)
img_raster = png::readPNG(Args$img)
img_grob <- grid::rasterGrob(img_raster, interpolate = FALSE, width = grid::unit(1,"npc"), height = grid::unit(1, "npc"))
object_spatial@misc$info[[`samplename`]]$img = img_64
object_spatial@misc$info[[`samplename`]]$img_gg = img_grob
object_spatial@misc$info[[`samplename`]]$size_x = as.integer(Args$size_x)
object_spatial@misc$info[[`samplename`]]$size_y = as.integer(Args$size_y)

mcsaveRDS(object_spatial, paste0(sampleid, ".rds"))
## step 3. plot


p <- FeaturePlot(object_spatial, 
    reduction = "spatial", pt.size = .05, 
    features = 'n_fragment', 
    min.cutoff = 'q05', max.cutoff = 'q95') +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(x = "", y = "")

ggsave(p, filename = paste0(sampleid, ".n_fragment.pdf"), 
    dpi = 600, width = 300, height = 120, units = "mm", limitsize = FALSE)

