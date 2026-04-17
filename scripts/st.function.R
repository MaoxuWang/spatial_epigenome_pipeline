library(patchwork)

extra_r_dir <- Sys.getenv("DSCHIC_R_FUNCTION_DIR", unset = "")
if (nzchar(extra_r_dir)) {
    source(file.path(extra_r_dir, "scRNA.func.R"))
    source(file.path(extra_r_dir, "plot_function.R"))
}


SpatialColors <- colorRampPalette(colors = c("black","#011481","#4933FB","#B932D9","yellow"), bias = 2)(10)

###color palettes
magenta_scale <- c("gray93", "#FFF7F3","#FDE0DD","#FCC5C0","#FA9FB5","#F768A1","#DD3497","#AE017E")
cyan_scale <- c("gray93","#F7FCF0","#E0F3DB","#CCEBC5","#7BCCC4","#4EB3D3","#2B8CBE","#0868AC","#084081","#084081","#084081")
fun_color_range <- colorRampPalette(c("gray93", "#AE017E"))  # Create color generating function
hm_colors <- fun_color_range(100)                         

#embryo
bright_visium <- c("#4477AA",  "#9ECAE1","#669C39", "#58BEC8", "#ccddaa", "#CCBB44", "#2F9558", "#B73D77","#E05B77", "#B06992", "#A8A8A8", "#E08762")
visium_dca <- c("#4477AA","#58BEC8","#669C39","#9ECAE1","#ccddaa","#CCBB44","#2F9558","#cc6677","#E05B77","#B73D77","#B06992","#A8A8A8","#E08762")
color_sections <- c("#E69F00", "#D55E00", "#aaaa00", "#009E73","#56B4E9", "#0072B2")
cols <-  c("#4477AA","#58AAD2","#58BEC8", "#2F9558","#669C39","#CCBB44", "#A8A8A8","#E05B77","#B73D77", "#B06992","#E08762") #embryo clusters res 0.7
cols_dca <-  c("#4477AA","#58AAD2","#58BEC8", "#2F9558","#669C39","#CCBB44", "#E05B77","#B73D77", "#B06992","#E08762")
cols_nmf <-c("#4477AA", "#66CCEE", "#228833", "#76A13B", "#CCBB44", "#BBBBBB","#EE6677", "#CB4C77", "#B27699","#AA3377", "#E69F00", "#D55E00")

#brca
colors_okate <- c("#E69F00",  "#F0E442", "#009E73",  "#56B4E9",  "#CC79A7")


background_vs_tissue <- function(
    all_obj,
    in_tissue_obj
){
    in_tissue_obj@meta.data %>%
        mutate(tissue = "tissue") %>%
        select(cellbarcode, tissue) -> in_tissue.df.some

    all_obj@meta.data %>%
        left_join(in_tissue.df.some) %>%
        mutate(tissue = if_else(is.na(tissue), "background", "tissue")) -> all.df

    plot <- ggplot(all.df, aes(x=log10(n_fragment), fill=tissue)) +
        geom_density(alpha = 0.8) + 
        scale_fill_manual(values = c("background" = "orange",
                                    "tissue" = "firebrick3"
                                    )) +
        theme_bw() +
        theme(  legend.title=element_blank(),
                axis.title.x =element_text(size=18), 
                axis.title.y =element_text(size=18),
                axis.text.x = element_text(size=10,colour = "black"),
                axis.text.y = element_text(size=10,colour = "black"),
                panel.grid = element_line(colour = 'white'))
    return(plot)
}

spatial_feature_plot <- function(
    seurat_obj,
    features,
    assay = "RNA",
    library_type = "DBiT", # visium-HD;Slide-tags
    pt.size.factor = 4,
    max_cutoff = 'q95',
    min_cutoff = 'q0',
    image.alpha = 0,
    cols = inferno(50),
    main_title = NULL,
    shape = 21,
    ...
){
    if (library_type == "DBiT"){
        stroke = 0.2
        shape = 22
    }else{
        stroke = 0
    }

    DefaultAssay(seurat_obj) <- assay

    if (length(features) > 1){

        plot_list <- lapply(seq_along(features), function(i) {  
            if (is.null(main_title)){
                main_title = features[i]
            }
            SpatialFeaturePlot(seurat_obj,
                features = features[i],
                image.alpha = image.alpha, 
                stroke = stroke, 
                max.cutoff = max_cutoff,
                min.cutoff = min_cutoff,
                pt.size.factor = pt.size.factor,
                shape = shape,
                ...) +
                theme(
                    legend.position = "right", 
                    legend.text=element_text(size=15), 
                    legend.title=element_text(size=15)) +
                scale_fill_gradientn(colours = cols) +
                fixed_coordinates +
                labs(title = main_title) +
                theme(
                            # legend.position="none",
                    plot.title = element_text(
                        hjust = 0.5, size = 15, face = "bold"),
                    legend.title = element_blank(), 
                    legend.key=element_rect(fill='white'), 
                    legend.text = element_text(size=12))
            })
    }else{
        if (is.null(main_title)){
            main_title = features
        }
        plot_list <- SpatialFeaturePlot(seurat_obj,
            features = features,
            image.alpha = image.alpha,
            stroke = stroke,
            max.cutoff = max_cutoff,
            min.cutoff = min_cutoff,
            pt.size.factor = pt.size.factor,
            shape = shape,
           ...) +
            theme(
                legend.position = "right",
                legend.text=element_text(size=15),
                legend.title=element_text(size=15)) +
            scale_fill_gradientn(colours = cols) +
            fixed_coordinates +
            labs(title = main_title) +
            theme(
                        # legend.position="none",
                plot.title = element_text(
                    hjust = 0.5, size = 15, face = "bold"),
                legend.title = element_blank(), 
                legend.key=element_rect(fill='white'), 
                legend.text = element_text(size=12))
    }
    return(plot_list)

}



spatial_dim_plot <- function(
    seurat_obj,
    column,
    group.by,
    library_type = "DBiT", # visium-HD;Slide-tags
    pt.size.factor = 3,
    image.alpha = 0,
    cols = NULL,
    pal = jdb_palette("corona"),
    main_title = NULL,
    shape = 21,
    label = F,
    ...
){
    if (library_type == "DBiT"){
        stroke = 0.2
        shape = 22
    }else{
        stroke = 0
    }
    
    if( ! is.factor(column)){
        column = factor(column, levels = unique(column))
        print("Please provide column information as factor")
    }
    if (is.null(cols)){
        # pal <- jdb_palette("lawhoops")
        if(length(unique(column)) > 30){
            cols = colorRampPalette(pal)(length(unique(column)))
        }else{
            cols <- pal[1:length(unique(column))]
        }
    }
    names(cols) <- levels(column)

    if (is.null(main_title)){
        main_title = group.by
    }

    Idents(seurat_obj) <- group.by
    p <- SpatialDimPlot(
        seurat_obj, 
        group.by = group.by, 
        label = label, repel = T, 
        pt.size.factor = pt.size.factor,
        image.alpha = image.alpha,
        shape = shape,
        stroke = stroke, 
        cols = cols,
        ...
        ) +
        theme(
            legend.position = "right", 
            legend.text=element_text(size=15), 
            legend.title=element_text(size=15)) +
        labs(title = main_title) +
        fixed_coordinates +
        theme(
                    # legend.position="none",
            plot.title = element_text(
                hjust = 0.5, size = 15, face = "bold"),
            legend.title = element_blank(), 
            # legend.key = element_rect(fill = "white", color = "white", linewidth = 2),
            legend.key.size = unit(8, "mm"),
            legend.text = element_text(size=12))

    return(p)

}  


# extractTFnames <-  function(motifIDs){
#     sapply(strsplit(sapply(strsplit(motifIDs,"_LINE.",fixed=FALSE),"[[",2),"_",fixed=FALSE),"[[",2)
#   }


addMotifAssay <- function(
    spatial_object,
    motif_data
){
    motif_deviation <- mcreadRDS(motif_data)
    ## extract tf names
    tf_names <- sapply(
        strsplit(rownames(motif_deviation), "_"), 
        function(x) x[3]
    )


    # rownames(motif_deviation) <- extractTFnames(rownames(motif_deviation))
    rownames(motif_deviation) <- tf_names

    motif_deviation.sub <- as(motif_deviation[,Cells(spatial_object)], "dgCMatrix")

    motif_assay <- CreateAssay5Object(
        data = motif_deviation.sub,  
        counts = NULL)
    spatial_object[['motif_activity']] <- motif_assay

    return(spatial_object)
}

spatial_clustering_banksy <- function(
    seurat_obj,
    assay = "RNA",
    resolution = 0.4,
    lambda = 0.6,
    k_geom = 15
){
    DefaultAssay(seurat_obj) <- assay

    seurat_obj <- RunBanksy(seurat_obj,
        lambda = 0.6, verbose = FALSE,
        assay = assay, slot = "data", features = "variable",
        k_geom = k_geom
    )

    DefaultAssay(seurat_obj) <- "BANKSY"
    seurat_obj %>%
        RunPCA(assay = "BANKSY", 
            reduction.name = "pca.banksy", 
            features = rownames(seurat_obj), 
            npcs = 30, verbose = F) %>%
        FindNeighbors(reduction = "pca.banksy", dims = 1:20, verbose = F) %>%
        FindClusters(cluster.name = "banksy_cluster", resolution = resolution, verbose = F) ->  seurat_obj

    return(seurat_obj)

}