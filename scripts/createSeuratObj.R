library(SeuratDisk)
library(Seurat)
library(tidyverse)
library(trqwe)


args = commandArgs(T)
input_h5ad = args[1] # .gene_activity.h5ad
metadata_file = args[2] # .metadata.txt
png_path = args[3] # cell_split/he_roi_small.png
barcode_pos_path = args[4] # 'super_spot/level_7/barcodes_pos.tsv.gz'
outdir = args[5]
radius_level = as.numeric(args[6])
in_tissue = as.logical(args[7]) # true or false

spot_radius = NULL
min.cells = 5
min.features = 100

sampleid = gsub(".gene_activity.h5ad", "", input_h5ad) # out/sample_name
sampleName = gsub(".*\\/", "", sampleid)
sampleid = paste0(outdir,"/", sampleName)


print(paste0("radius_level: ", radius_level))

#spot radius
spot_radius_lib <- c(0.00063, 0.00179, 0.0027, 0.0039, 0.004, 0.0045, 0.005, NA, NA, NA, NA, NA, 0.0120)
if(is.null(spot_radius)){
    spot_radius <- spot_radius_lib[radius_level]
}else{
    spot_radius = spot_radius
}

cal_zoom_rate <- function(width, height){
    std_width = 1000
    std_height = std_width / (46 * 31) * (46 * 36 * sqrt(3) / 2.0)
    if (std_width / std_height > width / height){
        scale = width / std_width
    }
    else{
        scale = height / std_height
    }
    return(scale)
}

#read barcode pos file
ReadBarcodePos <- function(barcode_pos_path){
    barcode_pos <- read.table(gzfile(barcode_pos_path),header = F) %>%
    dplyr::rename(Barcode = V1 , pos_w = V2, pos_h = V3)
    return(barcode_pos)
}


createSpatialObj <- function(
    spot_radius,
    object_spatial,
    png_path,
    barcode_pos_path){
    print(paste("spot_radius: ", spot_radius))
    ## add image 
    png <- png::readPNG(png_path)
    zoom_scale <-  cal_zoom_rate(dim(png)[2], dim(png)[1])

    ## get barcode pos file path
    barcode_pos <- ReadBarcodePos(barcode_pos_path = barcode_pos_path)
    barcode_pos <- barcode_pos %>% dplyr::filter(., Barcode %in% rownames(object_spatial@meta.data))

    ## make spatial coord file for seurat S4 class
    coord <- data.frame(tissue = 1,
                    row = barcode_pos$pos_h,
                    col = barcode_pos$pos_w,
                    imagerow = barcode_pos$pos_h,
                    imagecol = barcode_pos$pos_w)
    rownames(coord) <- barcode_pos$Barcode

    sample1 <-  new(Class = "VisiumV1",
                image = png,
                scale.factors = Seurat::scalefactors(zoom_scale, 100, zoom_scale, zoom_scale),
                coordinates = coord,
                spot.radius = spot_radius,
                assay = 'Spatial',
                key = "sample1_")
    object_spatial@images <- list(sample1 = sample1)

    object_spatial %>%
        SCTransform(assay = 'Spatial') -> object_spatial 
    
    object_spatial@images$sample1@spot.radius = spot_radius
    print(paste("spot_radius in seurat: ", object_spatial@images$sample1@spot.radius))

    return(object_spatial)
}


plot_spatial <- function(
    object_spatial,
    prefix){
    print(paste("spot_radius in seurat: ", object_spatial@images$sample1@spot.radius))

    QC_plot <- SpatialFeaturePlot(object_spatial,
        image.alpha = 0,
        features = c("nCount_Spatial",
        "nFeature_Spatial", 
        "n_fragment",
        "frac_dup",
        "tsse"), ncol = 5) & theme(legend.position = "right")

    marker_plot <- SpatialFeaturePlot(object_spatial,
                    image.alpha = 0,
                    features = c("Slc17a7", "Gad2", #  Glutamatergic neuron and GABAergic neuron     
                                'Arpp21', 'Nrg3', # Excitatory Neuron
                                    'Snhg11', 'Pvalb', # Inhibitory Neuron
                                    'Slc1a3', 'Slc1a2', # Astrocyte
                                    'Rgs9', 'Rarb', # Medium_spiny_neurons
                                    'Tgfbr1', 'Srgap2', # Microglial_cells
                                    'Lhfpl3', 'Dscam',  # Oligodendrocyte_precursor_cells
                                    'Mbp', 'Mog' # Oligodendrocyte
                                ),ncol = 2,
                    )& theme(legend.position = "right")

    ggsave(QC_plot, filename = paste0(prefix, ".QC.pdf"), dpi = 300, units = "mm", width = 400, height = 80, limitsize = FALSE)

    mcsaveRDS(object = object_spatial, file = paste0(prefix, ".spatial.rds"))
    print("rds file saved.")
    ggsave(marker_plot, filename = paste0(prefix, ".marker_gene.pdf"), dpi = 300, width = 200, height = 600, units = "mm", limitsize = FALSE)
}


## Step 1. convert h5ad to seurat

Convert(input_h5ad, 
        meta.data = FALSE,
        misc = FALSE,
        dest = "h5seurat", overwrite = TRUE)

input_h5seurat <- gsub("h5ad", "h5seurat", input_h5ad)

object <- LoadH5Seurat(input_h5seurat, meta.data = F)

metadata <- read.table(metadata_file, 
    sep = "\t", 
    row.names = 1, 
    header = TRUE)
metadata$cellbarcode = rownames(metadata)

assayData = as.matrix(object@assays$RNA@counts)
object_spatial <- Seurat::CreateSeuratObject(counts = assayData,
                           assay = 'Spatial',
                           min.cells=min.cells,
                           min.features=min.features)
object_spatial@meta.data %>%
    mutate(cellbarcode = rownames(object_spatial@meta.data)) %>%
    inner_join(metadata) -> object_spatial@meta.data 
rownames(object_spatial@meta.data) <- object_spatial@meta.data$cellbarcode

## Step 2. make spatial object
object_spatial.all <- createSpatialObj(
    spot_radius,
    object_spatial,
    png_path,
    barcode_pos_path)

plot_spatial(object_spatial.all, sampleid)

## in-tissue
if (in_tissue){
    in_tissue_path = gsub("barcodes_pos.tsv.gz", "barcodes_in_tissue.tsv.gz", barcode_pos_path)
    print(paste0("Extracting barcodes in tissue :", in_tissue_path))

    barcode_in_tissue = read_delim(in_tissue_path, delim = "\t", col_names = F)

    object_spatial.in_tissue <- subset(object_spatial, cells = barcode_in_tissue$X1)
    barcode_pos_path.in_tissue <- gsub("barcodes_pos.tsv.gz", "barcodes_in_tissue_pos.tsv.gz", barcode_pos_path)
    print(paste0("Position of barcodes in tissue from :", barcode_pos_path.in_tissue))

    object_spatial.in_tissue <- createSpatialObj(
        spot_radius,
        object_spatial.in_tissue,
        png_path,
        barcode_pos_path.in_tissue)

    plot_spatial(object_spatial.in_tissue, paste0(sampleid, "_in_tissue"))
}