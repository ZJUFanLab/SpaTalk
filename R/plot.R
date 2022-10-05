#' @title Plot spatial transcriptomics data
#'
#' @description Plot scatterpie for spot-based ST data
#' @param st_meta st_meta generated from \code{\link{generate_spot}}
#' @param pie_scale Scale of each pie to plot. Default is \code{1}.
#' @param xy_ratio Ratio of y and x coordinates. Default is \code{1}.
#' @param color Filled of colors for pie plot, length of \code{color} must be equal to the number of unique cell types in \code{sc_celltype}.
#' @import ggplot2 scatterpie
#' @importFrom ggpubr get_palette
#' @export

plot_st_pie_generate <- function(st_meta, pie_scale = 1, xy_ratio = 1, color = NULL) {
    # check
    if (is(st_meta, "SpaTalk")) {
        stop("This function is for st_meta generated from 'generate_spot'. Use 'generate_spot' instead to plot SpaTalk object!")
    }
    st_meta$x <- (st_meta$x - min(st_meta$x))/(max(st_meta$x) - min(st_meta$x))
    st_meta$y <- (st_meta$y - min(st_meta$y))/(max(st_meta$y) - min(st_meta$y))
    st_meta$y <- st_meta$y * xy_ratio
    cellname <- colnames(st_meta)[-c(1:4)]
    if (is.null(color)) {
        col_manual <- ggpubr::get_palette(palette = "lancet", k = length(cellname))
    } else {
        col_manual <- color
    }
    ggplot2::ggplot() + scatterpie::geom_scatterpie(data = st_meta, ggplot2::aes(x = x, y = y), col = cellname,
        color = NA, pie_scale = pie_scale) + ggplot2::coord_fixed(ratio = 1) + ggplot2::scale_fill_manual(values = col_manual) +
        ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank()) + ggplot2::labs(x = "scaled_x",
        y = "scaled_y")
}

#' @title Plot spatial transcriptomics data
#'
#' @description Plot scatterpie for spatial transcriptomics data
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param pie_scale Scale of each pie to plot. Default is \code{1}.
#' @param xy_ratio Ratio of y and x coordinates. Default is \code{1}.
#' @param color Filled of colors for pie plot, length of \code{color} must be equal to the number of unique cell types in \code{sc_celltype}.
#' @import ggplot2 scatterpie
#' @importFrom ggpubr get_palette
#' @export

plot_st_pie <- function(object, pie_scale = 1, xy_ratio = 1, color = NULL) {
    # check
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    st_type <- object@para$st_type
    if (st_type == "single-cell") {
        stop("This function is only for spot-based spatial transcriptomics data!")
    }
    st_meta <- object@meta$rawmeta
    st_meta$x <- (st_meta$x - min(st_meta$x))/(max(st_meta$x) - min(st_meta$x))
    st_meta$y <- (st_meta$y - min(st_meta$y))/(max(st_meta$y) - min(st_meta$y))
    st_meta$y <- st_meta$y * xy_ratio
    cellname <- colnames(st_meta)[-c(1:7)]
    if (is.null(color)) {
        col_manual <- ggpubr::get_palette(palette = "lancet", k = length(cellname))
    } else {
        col_manual <- color
    }
    ggplot2::ggplot() + scatterpie::geom_scatterpie(data = st_meta, ggplot2::aes(x = x, y = y), col = cellname,
        color = NA, pie_scale = pie_scale) + ggplot2::coord_fixed(ratio = 1) + ggplot2::scale_fill_manual(values = col_manual) +
        ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank()) + ggplot2::labs(x = "scaled_x",
        y = "scaled_y")
}

#' @title Plot spatial distribution of gene
#'
#' @description Point plot with spatial distribution of a gene for transcriptomics data
#' @details Please set \code{if_use_newmeta} as \code{FALSE} to plot the spatial distribution of gene before \code{\link{dec_celltype}} for spot-based data.
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param gene Symbol of gene, e.g., 'AKT1'.
#' @param size Point size. Default is \code{1}.
#' @param color_low Color for the lowest value.
#' @param color_mid Color for the middle value for using \code{scale_color_gradient2}. Default is \code{NULL}.
#' @param color_high Color for the highest value.
#' @param color_midpoint Value for the middle scale. Default is \code{NULL}.
#' @param if_use_newmeta Whether to use newmeta o plot the spatial distribution of gene after \code{\link{dec_celltype}} for spot-based data. Default is \code{TRUE}.
#' @param celltype gene in which celltype to plot. Default is \code{NULL}. Set \code{if_use_newmeta} TRUE when using this parameter.
#' @param if_plot_others Whether to plot other cells when to use defined \code{celltype}.
#' @import ggplot2
#' @importFrom stats median
#' @export

plot_st_gene <- function(object, gene, size = 1, color_low = "grey", color_mid = NULL, color_high = "blue",
    color_midpoint = NULL, if_use_newmeta = T, celltype = NULL, if_plot_others = T) {
    # check
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    # get result from dec_celltype
    st_type <- object@para$st_type
    if_skip_dec_celltype <- object@para$if_skip_dec_celltype
    if (st_type == "single-cell") {
        st_meta <- object@meta$rawmeta
        st_data <- object@data
        if (if_skip_dec_celltype) {
            st_data <- object@data$rawdata
        } else {
            st_data <- object@data$rawndata
        }
    }
    if (st_type == "spot") {
        if (if_use_newmeta) {
            if (if_skip_dec_celltype) {
                warning("if_use_newmeta is not used when skiping dec_celltype()!")
                st_meta <- object@meta$rawmeta
                st_data <- object@data$rawdata
            } else {
                st_meta <- object@meta$newmeta
                st_data <- object@data$newdata
            }
        } else {
            st_meta <- object@meta$rawmeta
            if (if_skip_dec_celltype) {
                st_data <- object@data$rawdata
            } else {
                st_data <- object@data$rawndata
            }
        }
    }
    if (length(gene) > 1) {
        stop("Please provide a single gene name each time!")
    }
    if (!gene %in% rownames(st_data)) {
        stop(paste0(gene, " is not in the st_data!"))
    }
    st_meta$gene <- as.numeric(st_data[gene, ])
    if (!is.null(celltype)) {
        if (!celltype %in% st_meta$celltype) {
            stop(paste0(celltype, " is not in st_meta!"))
        }
        if (if_use_newmeta) {
            if (if_plot_others) {
                st_meta[st_meta$celltype != celltype, ]$gene <- 0
            } else {
                st_meta <- st_meta[st_meta$celltype == celltype, ]
            }
        }
    }
    p <- ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(st_meta$x, st_meta$y, color = st_meta$gene),
        size = size) + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank()) + ggplot2::labs(x = "x",
        y = "y", color = gene)
    if (is.null(color_mid)) {
        p + ggplot2::scale_color_gradient(low = color_low, high = color_high)
    } else {
        if (is.null(color_midpoint)) {
            color_midpoint <- stats::median(st_meta$gene)
        }
        p + ggplot2::scale_color_gradient2(low = color_low, mid = color_mid, high = color_high, midpoint = color_midpoint)
    }
}

#' @title Plot spatial distribution of a single cell type
#'
#' @description Ponit plot with spatial distribution of a single predicted cell type for transcriptomics data
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param celltype Name of cell type in the \code{sc_celltype}.
#' @param size Point size. Default is \code{1}.
#' @param color_celltype Color for the celltype of interest.
#' @param color_others Color for the others.
#' @import ggplot2
#' @export

plot_st_celltype <- function(object, celltype, size = 1, color_celltype = "blue", color_others = "grey") {
    # check
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    # get result from dec_celltype
    st_type <- object@para$st_type
    if_skip_dec_celltype <- object@para$if_skip_dec_celltype
    if (st_type == "single-cell") {
        st_meta <- object@meta$rawmeta
    }
    if (st_type == "spot") {
        if (if_skip_dec_celltype) {
            st_meta <- object@meta$rawmeta
        } else {
            st_meta <- object@meta$newmeta
        }
    }
    if (length(celltype) > 1) {
        stop("Please provide a single cell type name each time!")
    }
    if (!celltype %in% st_meta$celltype) {
        stop(paste0(celltype, " is not in the st_meta!"))
    }
    st_meta[st_meta$celltype != celltype, ]$celltype <- "others"
    st_meta$celltype <- factor(st_meta$celltype, levels = c(celltype, "others"))
    ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(st_meta$x, st_meta$y, color = st_meta$celltype),
        size = size) + ggplot2::scale_color_manual(values = c(color_celltype, color_others)) + ggplot2::theme_bw() +
        ggplot2::theme(panel.grid = ggplot2::element_blank()) + ggplot2::labs(x = "x", y = "y", color = "celltype")
}

#' @title Plot spatial density of a single cell type
#'
#' @description Plot spatial density of a single predicted cell type for transcriptomics data
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param celltype Name of cell type in the \code{sc_celltype}.
#' @param type Select 'contour' or 'raster'.
#' @param if_plot_point Whether to plot points when type is 'contour'.
#' @param point_color Point color.
#' @param point_size Point size. Default is \code{1}.
#' @param color_low Color for the lowest value.
#' @param color_mid Color for the middle value for using \code{scale_color_gradient2}. Default is \code{NULL}.
#' @param color_high Color for the highest value.
#' @param color_midpoint Value for the middle scale. Default is \code{NULL}.
#' @param size Line size when type is 'contour'. Default is \code{1}.
#' @import ggplot2
#' @importFrom stats median
#' @export

plot_st_celltype_density <- function(object, celltype, type, if_plot_point = T, point_color = NULL, point_size = 1,
    color_low = "grey", color_mid = NULL, color_high = "blue", color_midpoint = NULL, size = 1) {
    # check
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    if (length(type) > 1 | !type %in% c("contour", "raster")) {
        stop("Please provide the correct type!")
    }
    # get result from dec_celltype
    st_type <- object@para$st_type
    if_skip_dec_celltype <- object@para$if_skip_dec_celltype
    if (st_type == "single-cell") {
        st_meta <- object@meta$rawmeta
    }
    if (st_type == "spot") {
        if (if_skip_dec_celltype) {
            st_meta <- object@meta$rawmeta
        } else {
            st_meta <- object@meta$newmeta
        }
    }
    if (length(celltype) > 1) {
        stop("Please provide a single cell type name each time!")
    }
    if (!celltype %in% st_meta$celltype) {
        stop(paste0(celltype, " is not in the st_meta!"))
    }
    if (is.null(point_color)) {
        cellname <- unique(st_meta$celltype)
        cellname <- cellname[order(cellname)]
        if ("unsure" %in% cellname) {
            cellname <- cellname[-which(cellname == "unsure")]
        }
        col_manual <- ggpubr::get_palette(palette = "lancet", k = length(cellname))
        point_color <- col_manual[which(cellname == celltype)]
    }
    st_meta <- st_meta[st_meta$celltype == celltype, ]
    p <- ggplot2::ggplot(data = st_meta, ggplot2::aes(x, y)) + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank())
    if (type == "contour") {
        if (if_plot_point) {
            p <- p + ggplot2::geom_point(size = point_size, color = point_color)
        }
        p <- p + ggplot2::stat_density2d(ggplot2::aes(color = ..level..), size = size)
        if (is.null(color_mid)) {
            p + scale_color_gradient(low = color_low, high = color_high)
        } else {
            if (is.null(color_midpoint)) {
                color_midpoint <- stats::median(st_meta$gene)
            }
            p + scale_color_gradient2(low = color_low, mid = color_mid, high = color_high, midpoint = color_midpoint)
        }
    } else {
        p <- p + ggplot2::stat_density2d(ggplot2::aes(fill = ..density..), geom = "raster", contour = F)
        if (is.null(color_mid)) {
            p + scale_fill_gradient(low = color_low, high = color_high)
        } else {
            if (is.null(color_midpoint)) {
                color_midpoint <- stats::median(st_meta$gene)
            }
            p + scale_fill_gradient2(low = color_low, mid = color_mid, high = color_high, midpoint = color_midpoint)
        }
    }
}

#' @title Plot spatial distribution of a single cell type percent
#'
#' @description Plot spatial distribution of a single predicted cell type percent for transcriptomics data
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param celltype Name of cell type in the \code{sc_celltype}.
#' @param size Point size. Default is \code{1}.
#' @param color_low Color for the lowest value.
#' @param color_mid Color for the middle value for using \code{scale_color_gradient2}. Default is \code{NULL}.
#' @param color_high Color for the highest value.
#' @param color_midpoint Value for the middle scale. Default is \code{NULL}.
#' @import ggplot2
#' @importFrom stats median
#' @export

plot_st_celltype_percent <- function(object, celltype, size = 1, color_low = NULL, color_mid = NULL, color_high = NULL,
    color_midpoint = NULL) {
    # check
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    if_skip_dec_celltype <- object@para$if_skip_dec_celltype
    if (if_skip_dec_celltype) {
        stop("Not availible when skiping dec_celltype()!")
    }
    # get result from dec_celltype
    st_meta <- object@meta$rawmeta
    if (length(celltype) > 1) {
        stop("Please provide a single cell type name each time!")
    }
    if (!celltype %in% colnames(st_meta)) {
        stop(paste0(celltype, " is not in the st_meta!"))
    }
    st_meta$celltype_percent <- as.numeric(st_meta[, celltype])
    if (is.null(color_low)) {
        color_low <- "grey"
    }
    if (is.null(color_high)) {
        color_high <- "blue"
    }
    if (is.null(color_mid)) {
        ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(st_meta$x, st_meta$y, color = st_meta$celltype_percent),
            size = size) + ggplot2::scale_color_gradient(low = color_low, high = color_high) + ggplot2::theme_bw() +
            ggplot2::theme(panel.grid = ggplot2::element_blank()) + ggplot2::labs(x = "x", y = "y", color = celltype)
    } else {
        if (is.null(color_midpoint)) {
            color_midpoint <- stats::median(st_meta$celltype_percent)
        }
        ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(st_meta$x, st_meta$y, color = st_meta$celltype_percent),
            size = size) + ggplot2::scale_color_gradient2(low = color_low, mid = color_mid, high = color_high,
            midpoint = color_midpoint) + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank()) +
            ggplot2::labs(x = "x", y = "y", color = paste0("Percent of ", celltype))
    }
}

#' @title Plot spatial distribution of all cell types
#'
#' @description Plot spatial distribution of all predicted cell types for transcriptomics data
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param size Point size. Default is \code{1}.
#' @param color Color for all predicted cell types.
#' @import ggplot2
#' @importFrom ggpubr get_palette
#' @export

plot_st_celltype_all <- function(object, size = 1, color = NULL) {
    # check
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    # get result from dec_celltype
    st_type <- object@para$st_type
    if_skip_dec_celltype <- object@para$if_skip_dec_celltype
    if (st_type == "single-cell") {
        st_meta <- object@meta$rawmeta
    }
    if (st_type == "spot") {
        if (if_skip_dec_celltype) {
            st_meta <- object@meta$rawmeta
        } else {
            st_meta <- object@meta$newmeta
        }
    }
    if (is.null(color)) {
        col_manual <- ggpubr::get_palette(palette = "lancet", k = length(unique(st_meta$celltype)))
    } else {
        col_manual <- color
    }
    ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(st_meta$x, st_meta$y, color = st_meta$celltype),
        size = size) + ggplot2::scale_color_manual(values = col_manual) + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank()) +
        ggplot2::labs(x = "x", y = "y", color = "celltype")
}

#' @title Plot heatpmap of correlation between marker genes and cell types
#'
#' @description Plot heatpmap of correlation between the expression of marker genes and the predicted score of cell types among all spatial cells or spots.
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param marker_genes A character containing the known marker genes to plot, provide at least two marker genes of interest.
#' @param celltypes A character containing name of cell type in the \code{sc_celltype}. Default is to plot all cell types.
#' @param color_low Color for the lowest value.
#' @param color_mid Color for the middle value for using \code{scale_color_gradient2}. Default is \code{NULL}.
#' @param color_high Color for the highest value.
#' @param scale Character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Corresponding values are 'row', 'column' and 'none'.
#' @param if_show_top Whether to plot a symbol to the highest value across rows or columns. Default is \code{TRUE}.
#' @param top_direction Direction to identify the highest value, select \code{'row'} or \code{'column'}.
#' @param border_color Color of the cell border. Default is \code{'NA'}.
#' @import pheatmap
#' @importFrom stats cor
#' @importFrom grDevices colorRampPalette
#' @importFrom ggpubr get_palette
#' @export

plot_st_cor_heatmap <- function(object, marker_genes, celltypes, color_low = NULL, color_mid = NULL, color_high = NULL,
    scale = "none", if_show_top = T, top_direction = "row", border_color = NA) {
    # check
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    if_skip_dec_celltype <- object@para$if_skip_dec_celltype
    if (if_skip_dec_celltype) {
        stop("Not availible when skiping dec_celltype()!")
    }
    # get result from dec_celltype
    st_meta <- object@meta$rawmeta
    st_data <- object@data$rawndata
    st_type_coef <- object@coef
    if (!all(celltypes %in% colnames(st_type_coef))) {
        celltypes <- celltypes[which(!celltypes %in% colnames(st_type_coef))]
        celltypes1 <- celltypes[1]
        if (length(celltypes) > 1) {
            for (i in 2:length(celltypes)) {
                celltypes1 <- paste(celltypes1, celltypes[i], sep = ", ")
            }
            stop(paste0(celltypes1, " are not in st_meta!"))
        } else {
            stop(paste0(celltypes1, " is not in st_meta!"))
        }
    }
    if (!all(marker_genes %in% rownames(st_data))) {
        marker_genes <- marker_genes[which(!marker_genes %in% rownames(st_data))]
        marker_genes1 <- marker_genes[1]
        if (length(marker_genes) > 1) {
            for (i in 2:length(marker_genes)) {
                marker_genes1 <- paste(marker_genes1, marker_genes[i], sep = ", ")
            }
            stop(paste0(marker_genes1, " are not in st_data!"))
        } else {
            stop(paste0(marker_genes1, " is not in st_data!"))
        }
    }
    st_data <- st_data[marker_genes, ]
    if (is.null(color_low)) {
        color_low <- "blue"
    }
    if (is.null(color_mid)) {
        color_mid <- "yellow"
    }
    if (is.null(color_high)) {
        color_high <- "red"
    }
    if (is.null(celltypes)) {
        celltypes <- colnames(st_type_coef)
    }
    plot_cor <- matrix(0, nrow = length(marker_genes), ncol = length(celltypes))
    rownames(plot_cor) <- marker_genes
    colnames(plot_cor) <- celltypes
    # colnames(plot_cor) <- celltypes
    for (i in 1:length(marker_genes)) {
        marker_genes1 <- marker_genes[i]
        marker_gene_data <- as.numeric(st_data[marker_genes1, ])
        for (j in 1:length(celltypes)) {
            celltypes1 <- celltypes[j]
            celltype_data <- as.numeric(st_meta[, celltypes1])
            plot_cor[i, j] <- stats::cor(marker_gene_data, celltype_data)
        }
    }
    heat_col <- (grDevices::colorRampPalette(c(color_low, color_mid, color_high)))(100)
    if (if_show_top) {
        if (top_direction == "row") {
            plot_cor1 <- plot_cor
            for (i in 1:nrow(plot_cor)) {
                data_cor <- as.numeric(plot_cor[i, ])
                plot_cor1[i, ] <- ""
                plot_cor1[i, which.max(data_cor)] <- "*"
            }
        } else {
            plot_cor1 <- plot_cor
            for (i in 1:ncol(plot_cor)) {
                data_cor <- as.numeric(plot_cor[, i])
                plot_cor1[, i] <- ""
                plot_cor1[which.max(data_cor), i] <- "*"
            }
        }
        pheatmap::pheatmap(mat = plot_cor, cluster_rows = F, cluster_cols = F, border_color = border_color,
            scale = scale, color = heat_col, display_numbers = plot_cor1)
    } else {
        pheatmap::pheatmap(mat = plot_cor, cluster_rows = F, cluster_cols = F, border_color = border_color,
            scale = scale, color = heat_col)
    }
}

#' @title Plot cell-cell distribution
#'
#' @description Point plot with spatial distribution of celltype_sender and celltype_receiver
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param celltype_sender Name of celltype_sender.
#' @param celltype_receiver Name of celltype_receiver.
#' @param color Color for celltype_sender, celltype_receiver, and others. Three values.
#' @param size Point size. Default is \code{1}.
#' @param if_plot_others Whether to plot others. Default is \code{TRUE}.
#' @param if_plot_density Whether to plot marginal density plots. Default is \code{TRUE}.
#' @param if_plot_edge Whether to plot edge between neighbors. Default is \code{TRUE}.
#' @param if_show_arrow Whether to show the arrow of the plotted edge. Default is \code{TRUE}.
#' @param arrow_length Arrow length.
#' @param plot_cells Which cells to plot. Default is all cells. Input a character vector of cell names to plot.
#' @import ggExtra ggplot2
#' @importFrom ggpubr get_palette
#' @importFrom ggExtra ggMarginal
#' @export

plot_ccdist <- function(object, celltype_sender, celltype_receiver, color = NULL, size = 1, if_plot_others = T,
    if_plot_density = T, if_plot_edge = T, if_show_arrow =T, arrow_length = 0.05, plot_cells = NULL) {
    # check
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    # get result from dec_celltype
    st_type <- object@para$st_type
    if_skip_dec_celltype <- object@para$if_skip_dec_celltype
    if (st_type == "single-cell") {
        st_meta <- object@meta$rawmeta
    } else {
        if (if_skip_dec_celltype) {
            st_meta <- object@meta$rawmeta
            colnames(st_meta)[1] <- "cell"
        } else {
            st_meta <- object@meta$newmeta
        }
    }
    if (!celltype_sender %in% st_meta$celltype) {
        stop("Please provide the correct name of celltype_sender!")
    }
    if (!celltype_receiver %in% st_meta$celltype) {
        stop("Please provide the correct name of celltype_receiver!")
    }
    if (!is.null(plot_cells)) {
        if (all(plot_cells %in% st_meta$cell)) {
            st_meta <- st_meta[st_meta$cell %in% plot_cells, ]
        } else {
            warning("Exclude the cells that is not in st_meta!")
            plot_cells <- plot_cells[plot_cells %in% st_meta$cell]
            if (length(plot_cells) < 10) {
                stop("Number of cells is less than 10!")
            }
            st_meta <- st_meta[st_meta$cell %in% plot_cells, ]
        }
    }
    cell_pair <- object@cellpair
    cell_pair <- cell_pair[[paste0(celltype_sender, " -- ", celltype_receiver)]]
    if (is.null(cell_pair)) {
        stop("No LR pairs found from the celltype_sender to celltype_receiver!")
    }
    cell_pair <- cell_pair[cell_pair$cell_sender %in% st_meta$cell & cell_pair$cell_receiver %in% st_meta$cell, ]
    cell_pair$x1 <- 0
    cell_pair$y1 <- 0
    cell_pair$x2 <- 0
    cell_pair$y2 <- 0
    for (i in 1:nrow(cell_pair)) {
        d1 <- cell_pair$cell_sender[i]
        d2 <- st_meta[st_meta$cell == d1, ]
        cell_pair$x1[i] <- d2$x
        cell_pair$y1[i] <- d2$y
        d1 <- cell_pair$cell_receiver[i]
        d2 <- st_meta[st_meta$cell == d1, ]
        cell_pair$x2[i] <- d2$x
        cell_pair$y2[i] <- d2$y
    }
    st_meta[!st_meta$celltype %in% c(celltype_sender, celltype_receiver), ]$celltype <- "Others"
    st_meta[st_meta$celltype == celltype_sender, ]$celltype <- paste0("Sender: ", celltype_sender)
    st_meta[st_meta$celltype == celltype_receiver, ]$celltype <- paste0("Receiver: ", celltype_receiver)
    if (if_plot_others) {
        st_meta$celltype <- factor(st_meta$celltype, levels = c(paste0("Sender: ", celltype_sender), paste0("Receiver: ", celltype_receiver), "Others"))
    } else {
        st_meta <- st_meta[st_meta$celltype != "Others", ]
        st_meta$celltype <- factor(st_meta$celltype, levels = c(paste0("Sender: ", celltype_sender), paste0("Receiver: ", celltype_receiver)))
    }
    if (is.null(color)) {
        cellname <- unique(as.character(st_meta$celltype))
        cellname <- cellname[order(cellname)]
        if ("unsure" %in% cellname) {
            cellname <- cellname[-which(cellname == "unsure")]
        }
        col_manual <- ggpubr::get_palette(palette = "lancet", k = length(cellname))
        if (if_plot_others) {
            color <- c(col_manual[which(cellname == paste0("Sender: ",celltype_sender))], col_manual[which(cellname == paste0("Receiver: ", celltype_receiver))], "grey80")
        } else {
            color <- c(col_manual[which(cellname == paste0("Sender: ",celltype_sender))], col_manual[which(cellname == paste0("Receiver: ", celltype_receiver))])
        }
    }
    p <- ggplot2::ggplot() + ggplot2::geom_point(data = st_meta, ggplot2::aes(x, y, color = celltype),
        size = size) + ggplot2::scale_color_manual(values = color) + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank())
    if (if_plot_edge) {
        if (if_show_arrow) {
            p <- p + geom_segment(data = cell_pair, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(arrow_length, "inches"), type = "closed"))
        } else {
            p <- p + geom_segment(data = cell_pair, aes(x = x1, y = y1, xend = x2, yend = y2))
        }
    }
    if (if_plot_density) {
        ggExtra::ggMarginal(p = p, type = "density", groupFill = T)
    } else {
        p
    }
}

#' @title Plot LR pairs
#'
#' @description Heatmap with LR pairs of celltype_sender and celltype_receiver
#' @param object SpaTalk object generated from \code{\link{dec_cci}}.
#' @param celltype_sender Name of celltype_sender.
#' @param celltype_receiver Name of celltype_receiver.
#' @param top_lrpairs Number of top lrpairs for plotting. Default is \code{20}.
#' @param color Color for the cells in heatmap.
#' @param border_color color of cell borders on heatmap, use NA if no border should be drawn.
#' @param type Set 'sig' to plot significant LRI pairs or set 'number' to plot the number of spatial LRI pairs.
#' @param fontsize_number fontsize of the numbers displayed in cells.
#' @param number_color color of the text.
#' @param color_low For 'number' type, define the color for the lowest value.
#' @param color_high For 'number' type, define the color for the highest value.
#' @import ggplot2 pheatmap grDevices
#' @importFrom ggpubr get_palette
#' @importFrom reshape2 dcast
#' @export

plot_cci_lrpairs <- function(object, celltype_sender, celltype_receiver, top_lrpairs = 20, color = NULL, border_color = "black", type = "sig",
    fontsize_number = 5, number_color = "black", color_low = NULL, color_high = NULL) {
    # check
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    # get result from dec_celltype
    st_type <- object@para$st_type
    if_skip_dec_celltype <- object@para$if_skip_dec_celltype
    if (st_type == "single-cell") {
        st_meta <- object@meta$rawmeta
        if (if_skip_dec_celltype) {
            st_data <- object@data$rawdata
        } else {
            st_data <- object@data$rawndata
        }
    } else {
        if (if_skip_dec_celltype) {
            st_meta <- object@meta$rawmeta
            colnames(st_meta)[1] <- "cell"
            st_data <- object@data$rawdata
        } else {
            st_meta <- object@meta$newmeta
            st_data <- object@data$newdata
        }
    }
    if (!celltype_sender %in% st_meta$celltype) {
        stop("Please provide the correct name of celltype_sender!")
    }
    if (!celltype_receiver %in% st_meta$celltype) {
        stop("Please provide the correct name of celltype_receiver!")
    }
    if (is.null(color)) {
        cellname <- unique(st_meta$celltype)
        cellname <- cellname[order(cellname)]
        if ("unsure" %in% cellname) {
            cellname <- cellname[-which(cellname == "unsure")]
        }
        col_manual <- ggpubr::get_palette(palette = "lancet", k = length(cellname))
        color <- col_manual[which(cellname == celltype_receiver)]
    }
    heat_col <- (grDevices::colorRampPalette(c("white", color)))(2)
    cell_pair <- object@cellpair
    cell_pair <- cell_pair[[paste0(celltype_sender, " -- ", celltype_receiver)]]
    if (is.null(cell_pair)) {
        stop("No LR pairs found from the celltype_sender to celltype_receiver!")
    }
    cell_pair <- cell_pair[cell_pair$cell_sender %in% st_meta$cell & cell_pair$cell_receiver %in% st_meta$cell, ]
    # get result from dec_cci
    lrpair <- object@lrpair
    lrpair <- lrpair[lrpair$celltype_sender == celltype_sender & lrpair$celltype_receiver == celltype_receiver, ]
    lrpair <- lrpair[order(-lrpair$score), ]
    if (nrow(lrpair) > top_lrpairs) {
        lrpair <- lrpair[1:top_lrpairs, ]
    }
    ligand <- unique(lrpair$ligand)
    receptor <- unique(lrpair$receptor)
    if (type == "sig") {
        # get lrpair_real
        lrpair_real <- object@lr_path$lrpairs
        lrpair_real <- lrpair_real[, c(1, 2)]
        lrpair_real$score <- 1
        lrpair_mat <- reshape2::dcast(lrpair_real, formula = ligand ~ receptor, fill = 0, value.var = "score")
        rownames(lrpair_mat) <- lrpair_mat$ligand
        lrpair_mat <- lrpair_mat[, -1]
        lrpair_mat <- lrpair_mat[ligand, receptor]
        if (!is.data.frame(lrpair_mat)) {
            stop("Limited number of ligand-receptor interactions!")
        }
        plot_res <- matrix("", nrow = length(ligand), ncol = length(receptor))
        rownames(plot_res) <- ligand
        colnames(plot_res) <- receptor
        for (i in 1:nrow(lrpair)) {
            plot_res[lrpair$ligand[i], lrpair$receptor[i]] <- "*"
        }
        pheatmap::pheatmap(lrpair_mat, cluster_cols = F, cluster_rows = F, color = heat_col, border_color = border_color, legend = F, display_numbers = plot_res,
            fontsize_number = fontsize_number, number_color = number_color, main = "Significantly enriched LRI pairs")
    } else {
        if (is.null(color_low)) {
            color_low <- "orange"
        }
        if (is.null(color_high)) {
            color_high <- "red"
        }
        lrpair_real <- lrpair[,c("ligand","receptor","lr_co_exp_num")]
        lrpair_mat <- reshape2::dcast(lrpair_real, formula = ligand ~ receptor, fill = 0, value.var = "lr_co_exp_num")
        rownames(lrpair_mat) <- lrpair_mat$ligand
        lrpair_mat <- lrpair_mat[, -1]
        heat_color <- grDevices::colorRampPalette(c(color_low, color_high))(max(as.matrix(lrpair_mat))-1)
        heat_color <- c("white", heat_color)
        pheatmap::pheatmap(lrpair_mat, cluster_cols = F, cluster_rows = F, border_color = border_color, color = heat_color, main = "Number of spatial LRI pairs")
    }
}

#' @title Plot LR pair
#'
#' @description Point plot with LR pair from celltype_sender to celltype_receiver
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param celltype_sender Name of celltype_sender.
#' @param celltype_receiver Name of celltype_receiver.
#' @param ligand Name of ligand from celltype_sender.
#' @param receptor Name of receptor from celltype_receiver.
#' @param color Color for ligand, receptor, and others. Three values.
#' @param size Point size. Default is \code{1}.
#' @param if_plot_density Whether to plot marginal density plots. Default is \code{TRUE}.
#' @param if_plot_edge Whether to plot edge between neighbors. Default is \code{TRUE}.
#' @param if_show_arrow Whether to show the arrow of the plotted edge. Default is \code{TRUE}.
#' @param arrow_length Arrow length.
#' @param plot_cells Which cells to plot. Default is all cells. Input a character vector of cell names to plot.
#' @import ggExtra ggplot2
#' @importFrom ggpubr get_palette
#' @importFrom ggExtra ggMarginal
#' @export

plot_lrpair <- function(object, celltype_sender, celltype_receiver, ligand, receptor, color = NULL, size = 1,
    if_plot_density = T, if_plot_edge = T, if_show_arrow =T, arrow_length = 0.05, plot_cells = NULL) {
    # check
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    # get result from dec_celltype
    st_type <- object@para$st_type
    if_skip_dec_celltype <- object@para$if_skip_dec_celltype
    if (st_type == "single-cell") {
        st_meta <- object@meta$rawmeta
        if (if_skip_dec_celltype) {
            st_data <- object@data$rawdata
        } else {
            st_data <- object@data$rawndata
        }
    } else {
        if (if_skip_dec_celltype) {
            st_meta <- object@meta$rawmeta
            colnames(st_meta)[1] <- "cell"
            st_data <- object@data$rawdata
        } else {
            st_meta <- object@meta$newmeta
            st_data <- object@data$newdata
        }
    }
    if (!celltype_sender %in% st_meta$celltype) {
        stop("Please provide the correct name of celltype_sender!")
    }
    if (!celltype_receiver %in% st_meta$celltype) {
        stop("Please provide the correct name of celltype_receiver!")
    }
    if (!ligand %in% rownames(st_data)) {
        stop("Please provide the correct name of ligand!")
    }
    if (!receptor %in% rownames(st_data)) {
        stop("Please provide the correct name of receptor!")
    }
    if (!is.null(plot_cells)) {
        if (all(plot_cells %in% st_meta$cell)) {
            st_meta <- st_meta[st_meta$cell %in% plot_cells, ]
        } else {
            warning("Exclude the cells that is not in st_meta!")
            plot_cells <- plot_cells[plot_cells %in% st_meta$cell]
            if (length(plot_cells) < 10) {
                stop("Number of cells is less than 10!")
            }
            st_meta <- st_meta[st_meta$cell %in% plot_cells, ]
        }
    }
    if (is.null(color)) {
        cellname <- unique(st_meta$celltype)
        cellname <- cellname[order(cellname)]
        if ("unsure" %in% cellname) {
            cellname <- cellname[-which(cellname == "unsure")]
        }
        col_manual <- ggpubr::get_palette(palette = "lancet", k = length(cellname))
        color <- c(col_manual[which(cellname == celltype_sender)], col_manual[which(cellname == celltype_receiver)], "grey80")
    }
    cell_pair <- object@cellpair
    cell_pair <- cell_pair[[paste0(celltype_sender, " -- ", celltype_receiver)]]
    if (is.null(cell_pair)) {
        stop("No LR pairs found from the celltype_sender to celltype_receiver!")
    }
    cell_pair <- cell_pair[cell_pair$cell_sender %in% st_meta$cell & cell_pair$cell_receiver %in% st_meta$cell, ]
    cell_pair$x1 <- 0
    cell_pair$y1 <- 0
    cell_pair$x2 <- 0
    cell_pair$y2 <- 0
    for (i in 1:nrow(cell_pair)) {
        d1 <- cell_pair$cell_sender[i]
        d2 <- st_meta[st_meta$cell == d1, ]
        cell_pair$x1[i] <- d2$x
        cell_pair$y1[i] <- d2$y
        d1 <- cell_pair$cell_receiver[i]
        d2 <- st_meta[st_meta$cell == d1, ]
        cell_pair$x2[i] <- d2$x
        cell_pair$y2[i] <- d2$y
    }
    st_data <- st_data[, st_meta$cell]
    st_meta$ligand <- as.numeric(st_data[ligand, ])
    st_meta$receptor <- as.numeric(st_data[receptor, ])
    st_meta$Expressed_genes <- "Others"
    st_meta_ligand <- st_meta[st_meta$celltype == celltype_sender & st_meta$ligand > 0, ]
    st_meta_receptor <- st_meta[st_meta$celltype == celltype_receiver & st_meta$receptor > 0, ]
    celltype_sender <- paste0(celltype_sender, ": ", ligand)
    celltype_receiver <- paste0(celltype_receiver, ": ", receptor)
    st_meta[st_meta$cell %in% st_meta_ligand$cell, ]$Expressed_genes <- celltype_sender
    st_meta[st_meta$cell %in% st_meta_receptor$cell, ]$Expressed_genes <- celltype_receiver
    cell_pair <- cell_pair[cell_pair$cell_sender %in% st_meta_ligand$cell, ]
    cell_pair <- cell_pair[cell_pair$cell_receiver %in% st_meta_receptor$cell, ]
    st_meta$Expressed_genes <- factor(st_meta$Expressed_genes, levels = c(celltype_sender, celltype_receiver, "Others"))
    p <- ggplot2::ggplot() + ggplot2::geom_point(data = st_meta, ggplot2::aes(x, y, color = Expressed_genes),
        size = size) + ggplot2::scale_color_manual(values = color) + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank())
    if (if_plot_edge) {
        if (if_show_arrow) {
            p <- p + geom_segment(data = cell_pair, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(arrow_length, "inches"), type = "closed"))
        } else {
            p <- p + geom_segment(data = cell_pair, aes(x = x1, y = y1, xend = x2, yend = y2))
        }
    }
    if (if_plot_density) {
        ggExtra::ggMarginal(p = p, type = "density", groupFill = T)
    } else {
        p
    }
}

#' @title Plot spatial distance of LR pair with vlnplot
#'
#' @description Violin plot spatial distance of LR pair between expressed senders and receivers and between expressed cell-cell pairs.
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param celltype_sender Name of celltype_sender.
#' @param celltype_receiver Name of celltype_receiver.
#' @param ligand Name of ligand from celltype_sender.
#' @param receptor Name of receptor from celltype_receiver.
#' @param vln_color Color for violins. Two values.
#' @param if_plot_boxplot Whether to plot boxplot. Default is \code{TRUE}.
#' @param box_width Box width. Default is \code{0.2}.
#' @import ggplot2
#' @importFrom ggpubr get_palette
#' @importFrom stats wilcox.test
#' @export

plot_lrpair_vln <- function(object, celltype_sender, celltype_receiver, ligand, receptor, vln_color = NULL,
    if_plot_boxplot = T, box_width = 0.2) {
    # check
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    # get result from dec_celltype
    st_type <- object@para$st_type
    if_skip_dec_celltype <- object@para$if_skip_dec_celltype
    if (st_type == "single-cell") {
        st_meta <- object@meta$rawmeta
        if (if_skip_dec_celltype) {
            st_data <- object@data$rawdata
        } else {
            st_data <- object@data$rawndata
        }
    } else {
        if (if_skip_dec_celltype) {
            st_meta <- object@meta$rawmeta
            colnames(st_meta)[1] <- "cell"
            st_data <- object@data$rawdata
        } else {
            st_meta <- object@meta$newmeta
            st_data <- object@data$newdata
        }
    }
    if (!celltype_sender %in% st_meta$celltype) {
        stop("Please provide the correct name of celltype_sender!")
    }
    if (!celltype_receiver %in% st_meta$celltype) {
        stop("Please provide the correct name of celltype_receiver!")
    }
    if (!ligand %in% rownames(st_data)) {
        stop("Please provide the correct name of ligand!")
    }
    if (!receptor %in% rownames(st_data)) {
        stop("Please provide the correct name of receptor!")
    }
    if (is.null(vln_color)) {
        cellname <- unique(st_meta$celltype)
        cellname <- cellname[order(cellname)]
        if ("unsure" %in% cellname) {
            cellname <- cellname[-which(cellname == "unsure")]
        }
        col_manual <- ggpubr::get_palette(palette = "lancet", k = length(cellname))
        vln_color <- c(col_manual[which(cellname == celltype_sender)], col_manual[which(cellname == celltype_receiver)], "grey80")
    }
    celltype_dist <- object@dist
    cell_pair <- object@cellpair
    cell_pair <- cell_pair[[paste0(celltype_sender, " -- ", celltype_receiver)]]
    if (is.null(cell_pair)) {
        stop("No LR pairs found from the celltype_sender to celltype_receiver!")
    }
    cell_pair <- cell_pair[cell_pair$cell_sender %in% st_meta$cell & cell_pair$cell_receiver %in% st_meta$cell, ]
    # calculate L-R distance across expressed cell-cell pairs
    ndata_ligand <- st_data[ligand, ]
    ndata_receptor <- st_data[receptor, ]
    ndata_ligand_cell <- names(ndata_ligand)[which(ndata_ligand > 0)]
    ndata_receptor_cell <- names(ndata_receptor)[which(ndata_receptor > 0)]
    celltype_dist1 <- celltype_dist[ndata_receptor_cell, ndata_ligand_cell]
    celltype_dist1 <- as.numeric(as.matrix(celltype_dist1))
    celltype_dist1 <- celltype_dist1[celltype_dist1 > 0]
    # calculate L-R distance across expressed senders and receivers
    ndata_ligand <- st_data[ligand, cell_pair$cell_sender]
    ndata_receptor <- st_data[receptor, cell_pair$cell_receiver]
    ndata_ligand_cell <- names(ndata_ligand)[which(ndata_ligand > 0)]
    ndata_receptor_cell <- names(ndata_receptor)[which(ndata_receptor > 0)]
    celltype_dist2 <- celltype_dist[ndata_receptor_cell, ndata_ligand_cell]
    celltype_dist2 <- as.numeric(as.matrix(celltype_dist2))
    celltype_dist2 <- celltype_dist2[celltype_dist2 > 0]
    if (length(celltype_dist1) > 0 & length(celltype_dist2) > 0) {
        print(stats::wilcox.test(celltype_dist1, celltype_dist2, alternative = "greater"))
    }
    d1 <- data.frame(celltype = "All cell-cell pairs", dist = celltype_dist1, stringsAsFactors = F)
    d2 <- data.frame(celltype = paste0(celltype_sender, "-", celltype_receiver, " pairs"), dist = celltype_dist2, stringsAsFactors = F)
    lr_dist_plot <- rbind(d1, d2)
    p <- ggplot(data = lr_dist_plot, aes(x = celltype, y = dist, fill = celltype)) + geom_violin(trim = T) + scale_fill_manual(values = vln_color)
    if (if_plot_boxplot) {
        p + geom_boxplot(width = box_width, outlier.shape = NA)
    }
}

#' @title Plot LR and downstream pathways
#'
#' @description Plot network with LR and downstream pathways
#' @param object SpaTalk object generated from \code{\link{dec_cci}}.
#' @param celltype_sender Name of celltype_sender.
#' @param celltype_receiver Name of celltype_receiver.
#' @param ligand Name of ligand from celltype_sender.
#' @param receptor Name of receptor from celltype_receiver.
#' @param color Color for points Two values.
#' @param size Size of points.
#' @param arrow_length Arrow length.
#' @import ggplot2 ggrepel
#' @importFrom ggpubr get_palette
#' @export

plot_lr_path <- function(object, celltype_sender, celltype_receiver, ligand, receptor, color = NULL, size = 5,
    arrow_length = 0.1) {
    # check
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    pathways <- object@lr_path$pathways
    ggi_tf <- unique(pathways[, c("src", "dest", "src_tf", "dest_tf")])
    # get result from dec_celltype
    st_type <- object@para$st_type
    if_skip_dec_celltype <- object@para$if_skip_dec_celltype
    if (st_type == "single-cell") {
        st_meta <- object@meta$rawmeta
        if (if_skip_dec_celltype) {
            st_data <- object@data$rawdata
        } else {
            st_data <- object@data$rawndata
        }
    } else {
        if (if_skip_dec_celltype) {
            st_meta <- object@meta$rawmeta
            colnames(st_meta)[1] <- "cell"
            st_data <- object@data$rawdata
        } else {
            st_meta <- object@meta$newmeta
            st_data <- object@data$newdata
        }
    }
    if (!celltype_sender %in% st_meta$celltype) {
        stop("Please provide the correct name of celltype_sender!")
    }
    if (!celltype_receiver %in% st_meta$celltype) {
        stop("Please provide the correct name of celltype_receiver!")
    }
    if (!ligand %in% rownames(st_data)) {
        stop("Please provide the correct name of ligand!")
    }
    if (!receptor %in% rownames(st_data)) {
        stop("Please provide the correct name of receptor!")
    }
    max_hop <- object@para$max_hop
    cell_pair <- object@cellpair
    cell_pair <- cell_pair[[paste0(celltype_sender, " -- ", celltype_receiver)]]
    if (is.null(cell_pair)) {
        stop("No LR pairs found from the celltype_sender to celltype_receiver!")
    }
    co_exp_ratio <- object@para$co_exp_ratio
    ggi_res <- .generate_ggi_res(ggi_tf, cell_pair, receptor, st_data, max_hop, co_exp_ratio)
    tf_gene_all <- .generate_tf_gene_all(ggi_res, max_hop)
    tf_gene_all <- data.frame(gene = names(tf_gene_all), hop = tf_gene_all, stringsAsFactors = F)
    tf_gene_all_new <- unique(tf_gene_all)
    tf_gene_all <- tf_gene_all_new$hop
    names(tf_gene_all) <- tf_gene_all_new$gene
    tf_path_all <- NULL
    for (i in 1:length(tf_gene_all)) {
        tf_path <- .get_tf_path(ggi_res, names(tf_gene_all)[i], as.numeric(tf_gene_all[i]), receptor)
        tf_path_all <- rbind(tf_path_all, tf_path)
    }
    tf_path_all$hop <- tf_path_all$hop + 2
    node_x <- unique(tf_path_all[, c("dest", "hop")])
    plot_node <- data.frame(gene = c(node_x$dest, ligand, receptor), x = c(node_x$hop, 1, 2), stringsAsFactors = F)
    node_y <- as.data.frame(table(plot_node$x), stringsAsFactors = F)
    node_y_max <- max(node_y$Freq)
    plot_node$y <- 0
    plot_node_new <- NULL
    for (i in 1:max(plot_node$x)) {
        plot_node_temp <- plot_node[plot_node$x == i, ]
        if (nrow(plot_node_temp) == node_y_max) {
            plot_node_temp$y <- 1:nrow(plot_node_temp)
        } else {
            node_y_inter <- (node_y_max - 1)/(nrow(plot_node_temp) + 1)
            node_y_new <- 1 + node_y_inter
            for (j in 1:nrow(plot_node_temp)) {
                plot_node_temp$y[j] <- node_y_new
                node_y_new <- node_y_new + node_y_inter
            }
        }
        plot_node_new <- rbind(plot_node_new, plot_node_temp)
    }
    plot_node_new$Expression <- 0
    plot_node_new$Celltype <- celltype_receiver
    st_data_gene <- st_data[plot_node_new$gene, st_meta[st_meta$celltype == celltype_receiver, ]$cell]
    plot_node_new$Expression <- as.numeric(rowMeans(as.matrix(st_data_gene)))
    st_data_ligand <- st_data[ligand, st_meta[st_meta$celltype == celltype_sender, ]$cell]
    plot_node_new[plot_node_new$gene == ligand, ]$Expression <- mean(st_data_ligand)
    plot_node_new[plot_node_new$gene == ligand, ]$Celltype <- celltype_sender
    plot_node_new$Celltype <- factor(plot_node_new$Celltype, levels = c(celltype_sender, celltype_receiver))
    if (is.null(color)) {
        cellname <- unique(st_meta$celltype)
        cellname <- cellname[order(cellname)]
        if ("unsure" %in% cellname) {
            cellname <- cellname[-which(cellname == "unsure")]
        }
        col_manual <- ggpubr::get_palette(palette = "lancet", k = length(cellname))
        col_manual <- c(col_manual[which(cellname == celltype_sender)], col_manual[which(cellname == celltype_receiver)])
    } else {
        col_manual <- color
    }
    tf_path_all <- tf_path_all[, c("src", "dest", "hop")]
    colnames(tf_path_all)[3] <- "dest_x"
    tf_path_all$src_x <- tf_path_all$dest_x - 1
    tf_path_all$src_y <- 0
    tf_path_all$dest_y <- 0
    for (i in 1:nrow(tf_path_all)) {
        gene <- tf_path_all$src[i]
        gene_x <- tf_path_all$src_x[i]
        plot_node_temp <- plot_node_new[plot_node_new$gene == gene & plot_node_new$x == gene_x, ]
        tf_path_all$src_y[i] <- plot_node_temp$y
        gene <- tf_path_all$dest[i]
        gene_x <- tf_path_all$dest_x[i]
        plot_node_temp <- plot_node_new[plot_node_new$gene == gene & plot_node_new$x == gene_x, ]
        tf_path_all$dest_y[i] <- plot_node_temp$y
    }
    tf_path_all <- tf_path_all[, c("src", "src_x", "src_y", "dest", "dest_x", "dest_y")]
    ligand_xy <- plot_node_new[plot_node_new$gene == ligand, ]
    ligand_xy <- ligand_xy[order(ligand_xy$x), ]
    ligand_x <- ligand_xy$x[1]
    ligand_y <- ligand_xy$y[1]
    receptor_xy <- plot_node_new[plot_node_new$gene == receptor, ]
    receptor_xy <- receptor_xy[order(receptor_xy$x), ]
    receptor_x <- receptor_xy$x[1]
    receptor_y <- receptor_xy$y[1]
    tf_path_temp <- data.frame(src = ligand, src_x = ligand_x, src_y = ligand_y, dest = receptor, dest_x = receptor_x, dest_y = receptor_y, stringsAsFactors = F)
    tf_path_all <- rbind(tf_path_all, tf_path_temp)
    plot_node_new$tf <- "NO"
    plot_node_new[plot_node_new$gene %in% names(tf_gene_all), ]$tf <- "YES"
    plot_node_new$Celltype <- factor(plot_node_new$Celltype, levels = c(celltype_sender, celltype_receiver))
    ggplot2::ggplot() + geom_segment(data = tf_path_all, aes(x = src_x, y = src_y, xend = dest_x, yend = dest_y)) +
        ggplot2::geom_point(data = plot_node_new, ggplot2::aes(x, y, color = Celltype), size = size) +
        ggplot2::scale_color_manual(values = col_manual) + ggrepel::geom_label_repel(data = plot_node_new,
        aes(x, y, label = gene, fill = tf)) + labs(x = NULL, y = NULL) + theme(axis.text = element_blank(),
        panel.grid = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white"))
}

#' @title River plot of significantly activated pathways and related downstream genes of receptors.
#'
#' @description River plot of significantly activated pathways and related downstream genes of receptors.
#' @param object SpaTalk object generated from \code{\link{dec_cci}}.
#' @param celltype_sender Name of celltype_sender.
#' @param celltype_receiver Name of celltype_receiver.
#' @param ligand Name of ligand from celltype_sender.
#' @param receptor Name of receptor from celltype_receiver.
#' @param min_gene_num Min genes number for each pathway.
#' @param pvalue P value of the Fisher-exact test.
#' @param color Color of pathways and genes. Two values.
#' @param color_flow Color of the flow.
#' @import ggalluvial ggplot2
#' @export

plot_path2gene <- function(object, celltype_sender, celltype_receiver, ligand, receptor, min_gene_num = 5,
    pvalue = 0.5, color = NULL, color_flow = "blue") {
    # check
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    tf_path_all <- get_lr_path(object, celltype_sender, celltype_receiver, ligand, receptor, min_gene_num)[[1]]
    pathways <- object@lr_path$pathways
    ggi_tf <- unique(pathways[, c("src", "dest", "src_tf", "dest_tf")])
    # get result from dec_celltype
    st_type <- object@para$st_type
    if_skip_dec_celltype <- object@para$if_skip_dec_celltype
    if (st_type == "single-cell") {
        st_meta <- object@meta$rawmeta
        if (if_skip_dec_celltype) {
            st_data <- object@data$rawdata
        } else {
            st_data <- object@data$rawndata
        }
    } else {
        if (if_skip_dec_celltype) {
            st_meta <- object@meta$rawmeta
            colnames(st_meta)[1] <- "cell"
            st_data <- object@data$rawdata
        } else {
            st_meta <- object@meta$newmeta
            st_data <- object@data$newdata
        }
    }
    if (!celltype_sender %in% st_meta$celltype) {
        stop("Please provide the correct name of celltype_sender!")
    }
    if (!celltype_receiver %in% st_meta$celltype) {
        stop("Please provide the correct name of celltype_receiver!")
    }
    if (!ligand %in% rownames(st_data)) {
        stop("Please provide the correct name of ligand!")
    }
    if (!receptor %in% rownames(st_data)) {
        stop("Please provide the correct name of receptor!")
    }
    # pathway
    ggi_pathway <- object@lr_path$pathways
    rec_pathway_all <- ggi_pathway[ggi_pathway$src == receptor | ggi_pathway$dest == receptor, ]
    if (nrow(rec_pathway_all) == 0) {
        stop("No significantly activated pathways found for this LRI")
    }
    rec_pathway_all <- unique(rec_pathway_all$pathway)
    rec_pathway_yes <- rep("NO", length(rec_pathway_all))
    rec_pathway_gene <- list()
    rec_pathway_pvalue <- as.double(rep(1, length(rec_pathway_all)))
    gene_rec <- unique(c(ligand, receptor, tf_path_all$src, tf_path_all$dest))
    gene_rec_num <- length(gene_rec)
    gene_all <- rownames(st_data)
    gene_all_num <- length(gene_all)
    for (j in 1:length(rec_pathway_all)) {
        gene_pathway <- ggi_pathway[ggi_pathway$pathway == rec_pathway_all[j], ]
        gene_pathway <- unique(c(gene_pathway$src, gene_pathway$dest))
        gene_pathway <- gene_pathway[gene_pathway %in% gene_all]
        if (length(gene_pathway) >= min_gene_num) {
            gene_pathway_num <- length(gene_pathway)
            gene_pathway_yes <- intersect(gene_rec, gene_pathway)
            gene_pathway_yes_num <- length(gene_pathway_yes)
            a <- matrix(c(gene_pathway_yes_num, gene_rec_num - gene_pathway_yes_num, gene_pathway_num -
                gene_pathway_yes_num, gene_all_num - gene_rec_num + gene_pathway_yes_num - gene_pathway_num), nrow = 2)
            fisher_pvalue <- fisher.test(a)
            fisher_pvalue <- as.double(fisher_pvalue$p.value)
            rec_pathway_gene[[j]] <- gene_pathway_yes
            rec_pathway_pvalue[j] <- fisher_pvalue
            rec_pathway_yes[j] <- "YES"
        }
    }
    gene2pathway_plot <- NULL
    for (j in 1:length(rec_pathway_all)) {
        if (rec_pathway_yes[j] == "YES" & rec_pathway_pvalue[j] < pvalue) {
            gene2pathway_plot_temp <- data.frame(pathway = rec_pathway_all[j], gene = rec_pathway_gene[[j]], stringsAsFactors = FALSE)
            gene2pathway_plot <- rbind(gene2pathway_plot, gene2pathway_plot_temp)
        }
    }
    if (is.null(gene2pathway_plot)) {
        stop("No significantly activated pathways found for this LRI")
    }
    gene2pathway_plot$num <- 1
    d1 <- gene2pathway_plot[, c("gene", "num")]
    d2 <- gene2pathway_plot[, c("pathway", "num")]
    d1$cluster <- "gene"
    d2$cluster <- "pathway"
    d1$group <- 1:nrow(d1)
    d2$group <- 1:nrow(d2)
    colnames(d1)[1] <- "type"
    colnames(d2)[1] <- "type"
    gene2pathway_plot <- rbind(d1, d2)
    gene_rec <- gene_rec[gene_rec %in% gene2pathway_plot$type]
    for (i in 1:length(gene_rec)) {
        d1 <- gene2pathway_plot[gene2pathway_plot$type == gene_rec[i], ]
        gene2pathway_plot[gene2pathway_plot$type == gene_rec[i], ]$num <- 1/nrow(d1)
    }
    rec_pathway_all <- rec_pathway_all[rec_pathway_all %in% gene2pathway_plot$type]
    coef <- length(gene_rec)/length(rec_pathway_all)
    for (i in 1:length(rec_pathway_all)) {
        d1 <- gene2pathway_plot[gene2pathway_plot$type == rec_pathway_all[i], ]
        gene2pathway_plot[gene2pathway_plot$type == rec_pathway_all[i], ]$num <- coef/nrow(d1)
    }
    p <- ggplot(gene2pathway_plot, aes(x = cluster, y = num, stratum = type, fill = cluster, alluvium = group,
        label = type)) + ggalluvial::geom_stratum(alpha = 1, width = 0.1) + ggalluvial::geom_flow(fill = color_flow,
        alpha = 0.25, width = 0.25) + labs(y = NULL) + geom_text(stat = "stratum", size = 3, angle = 0) +
        theme(axis.text.y.left = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(),
            panel.background = element_rect(fill = "white"))
    if (is.null(color)) {
        p
    } else {
        p + scale_fill_manual(values = color)
    }
}
