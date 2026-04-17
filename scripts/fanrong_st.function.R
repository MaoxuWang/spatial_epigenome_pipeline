
SpatialPlot_new <- function (object, group.by = NULL, features = NULL, images = NULL, 
          cols = NULL, image.alpha = 1, crop = TRUE, slot = "data", 
          min.cutoff = NA, max.cutoff = NA, cells.highlight = NULL, 
          cols.highlight = c("#DE2D26", "grey50"), facet.highlight = FALSE, 
          label = FALSE, label.size = 5, label.color = "white", 
          label.box = TRUE, repel = FALSE, ncol = NULL, combine = TRUE, 
          pt.size.factor = 1.6, alpha = c(1, 1), stroke = 0.25, interactive = FALSE, 
          do.identify = FALSE, identify.ident = NULL, do.hover = FALSE, 
          information = NULL) 
{
  if (isTRUE(x = do.hover) || isTRUE(x = do.identify)) {
    warning("'do.hover' and 'do.identify' are deprecated as we are removing plotly-based interactive graphics, use 'interactive' instead for Shiny-based interactivity", 
            call. = FALSE, immediate. = TRUE)
    interactive <- TRUE
  }
  if (!is.null(x = group.by) & !is.null(x = features)) {
    stop("Please specific either group.by or features, not both.")
  }
  images <- images %||% Images(object = object, assay = DefaultAssay(object = object))
  if (is.null(x = features)) {
    if (interactive) {
      return(ISpatialDimPlot(object = object, image = image, 
                             group.by = group.by, alpha = alpha))
    }
    group.by <- group.by %||% "ident"
    object[["ident"]] <- Idents(object = object)
    data <- object[[group.by]]
    for (group in group.by) {
      if (!is.factor(x = data[, group])) {
        data[, group] <- factor(x = data[, group])
      }
    }
  }
  else {
    if (interactive) {
      return(ISpatialFeaturePlot(object = object, feature = features[1], 
                                 image = images[1], slot = slot, alpha = alpha))
    }
    data <- FetchData(object = object, vars = features, slot = slot)
    features <- colnames(x = data)
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
      return(ifelse(test = is.na(x = cutoff), yes = min(data[, 
                                                             feature]), no = cutoff))
    }, cutoff = min.cutoff, feature = features)
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
      return(ifelse(test = is.na(x = cutoff), yes = max(data[, 
                                                             feature]), no = cutoff))
    }, cutoff = max.cutoff, feature = features)
    check.lengths <- unique(x = vapply(X = list(features, 
                                                min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check.lengths) != 1) {
      stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    data <- sapply(X = 1:ncol(x = data), FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <- SetQuantile(cutoff = min.cutoff[index], 
                             data.feature)
      max.use <- SetQuantile(cutoff = max.cutoff[index], 
                             data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      return(data.feature)
    })
    colnames(x = data) <- features
    rownames(x = data) <- Cells(x = object)
  }
  if (length(x = images) == 0) {
    images <- Images(object = object)
  }
  if (length(x = images) < 1) {
    stop("Could not find any spatial image information")
  }
  features <- colnames(x = data)
  colnames(x = data) <- features
  rownames(x = data) <- colnames(x = object)
  facet.highlight <- facet.highlight && (!is.null(x = cells.highlight) && 
                                           is.list(x = cells.highlight))
  if (do.hover) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning("'do.hover' requires only one image, using image ", 
              images, call. = FALSE, immediate. = TRUE)
    }
    if (length(x = features) > 1) {
      features <- features[1]
      type <- ifelse(test = is.null(x = group.by), yes = "feature", 
                     no = "grouping")
      warning("'do.hover' requires only one ", type, 
              ", using ", features, call. = FALSE, immediate. = TRUE)
    }
    if (facet.highlight) {
      warning("'do.hover' requires no faceting highlighted cells", 
              call. = FALSE, immediate. = TRUE)
      facet.highlight <- FALSE
    }
  }
  if (facet.highlight) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning("Faceting the highlight only works with a single image, using image ", 
              images, call. = FALSE, immediate. = TRUE)
    }
    ncols <- length(x = cells.highlight)
  }
  else {
    ncols <- length(x = images)
  }
  plots <- vector(mode = "list", length = length(x = features) * 
                    ncols)
  for (i in 1:ncols) {
    plot.idx <- i
    image.idx <- ifelse(test = facet.highlight, yes = 1, 
                        no = i)
    image.use <- object[[images[[image.idx]]]]
    coordinates <- GetTissueCoordinates(object = image.use)
    highlight.use <- if (facet.highlight) {
      cells.highlight[i]
    }
    else {
      cells.highlight
    }
    for (j in 1:length(x = features)) {
      cols.unset <- is.factor(x = data[, features[j]]) && 
        is.null(x = cols)
      if (cols.unset) {
        cols <- hue_pal()(n = length(x = levels(x = data[, 
                                                         features[j]])))
        names(x = cols) <- levels(x = data[, features[j]])
      }
      plot <- SingleSpatialPlot(data = cbind(coordinates, 
                                             data[rownames(x = coordinates), features[j], 
                                                  drop = FALSE]), image = image.use, image.alpha = image.alpha, 
                                col.by = features[j], cols = cols, alpha.by = if (is.null(x = group.by)) {
                                  features[j]
                                }
                                else {
                                  NULL
                                }, geom = if (inherits(x = image.use, what = "STARmap")) {
                                  "poly"
                                }
                                else {
                                  "spatial"
                                }, cells.highlight = highlight.use, cols.highlight = cols.highlight, 
                                pt.size.factor = pt.size.factor, stroke = stroke, 
                                crop = crop)
      if (is.null(x = group.by)) {
        
        plot <- plot + scale_fill_gradientn(name = features[j], 
                                            colours = c("blue","green", "red"), oob = scales::squish) + theme(legend.position = "top") + 
          scale_alpha(range = alpha) + guides(alpha = FALSE)
      }
      else if (label) {
        plot <- LabelClusters(plot = plot, id = ifelse(test = is.null(x = cells.highlight), 
                                                       yes = features[j], no = "highlight"), 
                              geom = if (inherits(x = image.use, what = "STARmap")) {
                                "GeomPolygon"
                              }
                              else {
                                "GeomSpatial"
                              }, repel = repel, size = label.size, color = label.color, 
                              box = label.box, position = "nearest")
      }
      if (j == 1 && length(x = images) > 1 && !facet.highlight) {
        plot <- plot + ggtitle(label = images[[image.idx]]) + 
          theme(plot.title = element_text(hjust = 0.5))
      }
      if (facet.highlight) {
        plot <- plot + ggtitle(label = names(x = cells.highlight)[i]) + 
          theme(plot.title = element_text(hjust = 0.5)) + 
          NoLegend()
      }
      plots[[plot.idx]] <- plot
      plot.idx <- plot.idx + ncols
      if (cols.unset) {
        cols <- NULL
      }
    }
  }
  if (length(x = images) > 1 && combine) {
    plots <- wrap_plots(plots = plots, ncol = length(x = images))
  }
  else if (length(x = images == 1) && combine) {
    plots <- wrap_plots(plots = plots, ncol = ncol)
  }
  return(plots)
}

SpatialDimPlot_new <- function (object, group.by = NULL, features = NULL, images = NULL, 
                             cols = NULL, image.alpha = 1, crop = TRUE, slot = "data", 
                             min.cutoff = NA, max.cutoff = NA, cells.highlight = NULL, 
                             cols.highlight = c("green", "grey50"), facet.highlight = FALSE, 
                             label = FALSE, label.size = 5, label.color = "white", 
                             label.box = TRUE, repel = FALSE, ncol = NULL, combine = TRUE, 
                             pt.size.factor = 1.6, alpha = c(1, 1), stroke = 0.25, interactive = FALSE, 
                             do.identify = FALSE, identify.ident = NULL, do.hover = FALSE, 
                             information = NULL) 
{
  if (isTRUE(x = do.hover) || isTRUE(x = do.identify)) {
    warning("'do.hover' and 'do.identify' are deprecated as we are removing plotly-based interactive graphics, use 'interactive' instead for Shiny-based interactivity", 
            call. = FALSE, immediate. = TRUE)
    interactive <- TRUE
  }
  if (!is.null(x = group.by) & !is.null(x = features)) {
    stop("Please specific either group.by or features, not both.")
  }
  images <- images %||% Images(object = object, assay = DefaultAssay(object = object))
  if (is.null(x = features)) {
    if (interactive) {
      return(ISpatialDimPlot(object = object, image = image, 
                             group.by = group.by, alpha = alpha))
    }
    group.by <- group.by %||% "ident"
    object[["ident"]] <- Idents(object = object)
    data <- object[[group.by]]
    for (group in group.by) {
      if (!is.factor(x = data[, group])) {
        data[, group] <- factor(x = data[, group])
      }
    }
  }
  else {
    if (interactive) {
      return(ISpatialFeaturePlot(object = object, feature = features[1], 
                                 image = images[1], slot = slot, alpha = alpha))
    }
    data <- FetchData(object = object, vars = features, slot = slot)
    features <- colnames(x = data)
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
      return(ifelse(test = is.na(x = cutoff), yes = min(data[, 
                                                             feature]), no = cutoff))
    }, cutoff = min.cutoff, feature = features)
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
      return(ifelse(test = is.na(x = cutoff), yes = max(data[, 
                                                             feature]), no = cutoff))
    }, cutoff = max.cutoff, feature = features)
    check.lengths <- unique(x = vapply(X = list(features, 
                                                min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check.lengths) != 1) {
      stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    data <- sapply(X = 1:ncol(x = data), FUN = function(index) {
      data.feature <- as.vector(x = data[, index])
      min.use <- SetQuantile(cutoff = min.cutoff[index], 
                             data.feature)
      max.use <- SetQuantile(cutoff = max.cutoff[index], 
                             data.feature)
      data.feature[data.feature < min.use] <- min.use
      data.feature[data.feature > max.use] <- max.use
      return(data.feature)
    })
    colnames(x = data) <- features
    rownames(x = data) <- Cells(x = object)
  }
  if (length(x = images) == 0) {
    images <- Images(object = object)
  }
  if (length(x = images) < 1) {
    stop("Could not find any spatial image information")
  }
  features <- colnames(x = data)
  colnames(x = data) <- features
  rownames(x = data) <- colnames(x = object)
  facet.highlight <- facet.highlight && (!is.null(x = cells.highlight) && 
                                           is.list(x = cells.highlight))
  if (do.hover) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning("'do.hover' requires only one image, using image ", 
              images, call. = FALSE, immediate. = TRUE)
    }
    if (length(x = features) > 1) {
      features <- features[1]
      type <- ifelse(test = is.null(x = group.by), yes = "feature", 
                     no = "grouping")
      warning("'do.hover' requires only one ", type, 
              ", using ", features, call. = FALSE, immediate. = TRUE)
    }
    if (facet.highlight) {
      warning("'do.hover' requires no faceting highlighted cells", 
              call. = FALSE, immediate. = TRUE)
      facet.highlight <- FALSE
    }
  }
  if (facet.highlight) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning("Faceting the highlight only works with a single image, using image ", 
              images, call. = FALSE, immediate. = TRUE)
    }
    ncols <- length(x = cells.highlight)
  }
  else {
    ncols <- length(x = images)
  }
  plots <- vector(mode = "list", length = length(x = features) * 
                    ncols)
  for (i in 1:ncols) {
    plot.idx <- i
    image.idx <- ifelse(test = facet.highlight, yes = 1, 
                        no = i)
    image.use <- object[[images[[image.idx]]]]
    coordinates <- GetTissueCoordinates(object = image.use)
    highlight.use <- if (facet.highlight) {
      cells.highlight[i]
    }
    else {
      cells.highlight
    }
    for (j in 1:length(x = features)) {
      cols.unset <- is.factor(x = data[, features[j]]) && 
        is.null(x = cols)
      if (cols.unset) {
        cols <- scales::hue_pal()(n = length(x = levels(x = data[, 
                                                         features[j]])))
        names(x = cols) <- levels(x = data[, features[j]])
      }
      plot <- SingleSpatialPlot(data = cbind(coordinates, 
                                             data[rownames(x = coordinates), features[j], 
                                                  drop = FALSE]), image = image.use, image.alpha = image.alpha, 
                                col.by = features[j], cols = cols, alpha.by = if (is.null(x = group.by)) {
                                  features[j]
                                }
                                else {
                                  NULL
                                }, geom = if (inherits(x = image.use, what = "STARmap")) {
                                  "poly"
                                }
                                else {
                                  "spatial"
                                }, cells.highlight = highlight.use, cols.highlight = cols.highlight, 
                                pt.size.factor = pt.size.factor, stroke = stroke, 
                                crop = crop)
      if (is.null(x = group.by)) {
        
        plot <- plot + scale_fill_gradientn(name = features[j], 
                                            colours = c("blue","green", "red"), oob = scales::squish) + theme(legend.position = "top") + 
          scale_alpha(range = alpha) + guides(alpha = FALSE)
      }
      else if (label) {
        plot <- LabelClusters(plot = plot, id = ifelse(test = is.null(x = cells.highlight), 
                                                       yes = features[j], no = "highlight"), 
                              geom = if (inherits(x = image.use, what = "STARmap")) {
                                "GeomPolygon"
                              }
                              else {
                                "GeomSpatial"
                              }, repel = repel, size = label.size, color = label.color, 
                              box = label.box, position = "nearest")
      }
      else if (!is.null(x = cells.highlight)){                 # add alpha for cell.highlight
        plot <- plot + scale_alpha_discrete(range = alpha)
      }
      
      if (j == 1 && length(x = images) > 1 && !facet.highlight) {
        plot <- plot + ggtitle(label = images[[image.idx]]) + 
          theme(plot.title = element_text(hjust = 0.5))
      }
      if (facet.highlight) {
        plot <- plot + ggtitle(label = names(x = cells.highlight)[i]) + 
          theme(plot.title = element_text(hjust = 0.5)) + 
          NoLegend()
      }
      plots[[plot.idx]] <- plot
      plot.idx <- plot.idx + ncols
      if (cols.unset) {
        cols <- NULL
      }
    }
  }
  if (length(x = images) > 1 && combine) {
    plots <- wrap_plots(plots = plots, ncol = length(x = images))
  }
  else if (length(x = images == 1) && combine) {
    plots <- wrap_plots(plots = plots, ncol = ncol)
  }
  return(plots)
}

