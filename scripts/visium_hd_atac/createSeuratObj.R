library(SeuratDisk)
library(Seurat)
library(tidyverse)
library(trqwe)
library(arrow)
library(argparse)
library(ggrastr)

parser <- ArgumentParser()
parser$add_argument("--input_h5ad", type="character", help="input h5ad file")
parser$add_argument("--metadata_file", type="character", help="metadata file")
parser$add_argument("--img_dir", type="character", help="image directory")
parser$add_argument("--outdir", type="character", help="output directory")
Args <- parser$parse_args()

spot_radius = NULL
min.cells = 5
min.features = 100

sampleid = gsub(".gene_activity.h5ad", "", Args$input_h5ad) # out/sample_name
sampleName = gsub(".*\\/", "", sampleid)
outdir = paste0(Args$outdir,"/", sampleName)

readImage <- function(
    img_dir){
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
        assay = "RNA",
        image.type = "VisiumV2",
        filter.matrix = TRUE)

    return(image)
}

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

object_spatial <- AddMetaData(object_spatial, metadata = metadata)

## Step 2. make spatial object

image <- readImage(Args$img_dir)

DefaultAssay(object = image) <- "RNA"

cell_used <- intersect(Cells(image), Cells(object_spatial))
image <- image[cell_used]
object_spatial <- subset(object_spatial, cells = cell_used)

object_spatial[['slice_1']] <- image
mcsaveRDS(object_spatial, paste0(outdir, ".rds"))
## step 3. plot

p <- SpatialFeaturePlot(object_spatial,
    features = "n_fragment",
    pt.size.factor = .1,
    stroke = 0,
    max.cutoff = "q98",
    min.cutoff = "q2",
    image.alpha = 0) 

p.rasterise <- ggrastr::rasterise(p, layers='Spatial',dpi = 300)

ggsave(p.rasterise, filename = paste0(outdir, ".n_fragment.pdf"), 
    dpi = 600, width = 120, height = 120, units = "mm", limitsize = FALSE)

