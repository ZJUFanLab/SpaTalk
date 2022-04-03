#' @title Definition of 'SpaTalk' class
#'
#' @description An S4 class containing the data, meta, and results of inferred cell type compositions, LR pairs, and pathways.
#' @slot data A list containing the raw and normalized data.
#' @slot meta A list containing the raw and new meta data.
#' @slot para A list containing the parameters.
#' @slot coef A matrix containing the results of deconvolution.
#' @slot cellpair A list containing the cell-cell pairs based on the spatial distance.
#' @slot dist A matrix containing the Euclidean distance among cells.
#' @slot lrpair A data frame containing the inferred LR pairs.
#' @slot tf A data frame containing the TFs of receptors.
#' @slot lr_path A list containing the lrpairs and pathways.
#' @import methods
#' @name SpaTalk
#' @rdname SpaTalk
#' @aliases SpaTalk-class
#' @exportClass SpaTalk

setClass("SpaTalk", representation(data = "list", meta = "list", para = "list", coef = "matrix",
    cellpair = "list", dist = "matrix", lrpair = "data.frame", tf = "data.frame", lr_path = "list"),
    prototype(data = list(), meta = list(), para = list(), coef = matrix(), cellpair = list(),
        dist = matrix(), lrpair = data.frame(), tf = data.frame(), lr_path = list()))

#' @title Decomposing cell type for spatial transcriptomics data
#'
#' @description Identify the cellular composition for single-cell or spot-based spatial transcriptomics data with non-negative regression.
#' @param object SpaTalk object generated from createSpaTalk.
#' @param sc_data A matrix containing counts of single-cell RNA-seq data as the reference, each column representing a cell, each row representing a gene.
#' @param sc_celltype A character containing the cell type of the reference single-cell RNA-seq data.
#' @param min_percent Min percent to predict new cell type for single-cell st_data or predict new cell for spot-based st_data. Default is \code{0.5}.
#' @param min_nFeatures Min number of expressed features/genes for each spot/cell in \code{st_data}. Default is \code{10}.
#' @param if_use_normalize_data Whether to use normalized \code{st_data} and \code{sc_data} with Seurat normalization. Default is \code{TRUE}.
#' @param if_use_hvg Whether to use highly variable genes for non-negative regression. Default is \code{FALSE}
#' @param if_use_all_cores Whether to use all CPU cores. Default is \code{TRUE}.
#' @param iter_num Number of iteration to genenrate the single-cell data for spot-based data. Default is \code{1000}.
#' @return SpaTalk object containing the decomposing results.
#' @import methods
#' @export

setGeneric("dec_celltype", def = function(object, sc_data, sc_celltype, min_percent = 0.5, min_nFeatures = 10,
    if_use_normalize_data = T, if_use_hvg = F, if_use_all_cores = T, iter_num = 1000) {
    standardGeneric("dec_celltype")
})

#' @title Find lrpairs and pathways
#'
#' @description Find \code{lrpairs} and \code{pathways} with receptors having downstream targets and transcriptional factors.
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param lrpairs A data.frame of the system data containing ligand-receptor pairs of \code{'Human'} and \code{'Mouse'} from CellTalkDB.
#' @param pathways A data.frame of the system data containing gene-gene interactions and pathways from KEGG and Reactome as well as the information of transcriptional factors.
#' @param max_hop Max hop from the receptor to the downstream target transcriptional factor to find for receiving cells. Default is \code{5}.
#' @return SpaTalk object containing the filtered lrpairs and pathways.
#' @import methods
#' @export

setGeneric("find_lr_path", def = function(object, lrpairs, pathways, max_hop = 5) {
    standardGeneric("find_lr_path")
})

#' @title Decomposing cell-cell communications for spatial transciptomics data
#'
#' @description Identify the cell-cell communications for single-cell or spot-based spatial transciptomics data with proximal ligand-receptor-target interactions.
#' @param object SpaTalk object after \code{\link{find_lr_path}}.
#' @param celltype_sender Name of celltype_sender.
#' @param celltype_receiver Name of celltype_receiver.
#' @param n_neighbor Number of neighbor cells to select as the proximal cell-cell pair. Default is \code{10}.
#' @param min_pairs Min proximal cell-cell pairs between for sending and receiving cell types. Default is \code{5}.
#' @param min_pairs_ratio Min proximal cell-cell pairs ratio between for sending and receiving cell types. Default is \code{0}.
#' @param per_num Number of repeat times for permutation test. Default is \code{1000}.
#' @param pvalue Include the significantly proximal LR pairs with this cutoff of p value from permutation test. Default is \code{0.05}.
#' @param co_exp_ratio Min cell ratio in receiving cells with co-expressed source and target genes for predicting the downstream pathway activity.
#' @return SpaTalk object containing the inferred LR pairs and pathways.
#' @import methods
#' @export

setGeneric("dec_cci", def = function(object, celltype_sender, celltype_receiver, n_neighbor = 10,
    min_pairs = 5, min_pairs_ratio = 0, per_num = 1000, pvalue = 0.05, co_exp_ratio = 0.1) {
    standardGeneric("dec_cci")
})

#' @title Decomposing cell-cell communications for spatial transciptomics data
#'
#' @description Identify the all cell-cell communications for single-cell or spot-based spatial transciptomics data with proximal ligand-receptor-target interactions.
#' @param object SpaTalk object after \code{\link{find_lr_path}}.
#' @param n_neighbor Number of neighbor cells to select as the proximal cell-cell pair. Default is \code{10}.
#' @param min_pairs Min proximal cell-cell pairs between for sending and receiving cell types. Default is \code{5}.
#' @param min_pairs_ratio Min proximal cell-cell pairs ratio between for sending and receiving cell types. Default is \code{0}.
#' @param per_num Number of repeat times for permutation test. Default is \code{1000}.
#' @param pvalue Include the significantly proximal LR pairs with this cutoff of p value from permutation test. Default is \code{0.05}.
#' @param co_exp_ratio Min cell ratio in receiving cells with co-expressed source and target genes for predicting the downstream pathway activity.
#' @return SpaTalk object containing the inferred LR pairs and pathways.
#' @import methods
#' @export

setGeneric("dec_cci_all", def = function(object, n_neighbor = 10, min_pairs = 5, min_pairs_ratio = 0,
    per_num = 1000, pvalue = 0.05, co_exp_ratio = 0.1) {
    standardGeneric("dec_cci_all")
})

#' @title Plot spatial transcriptomics data
#'
#' @description Plot scatterpie for spatial transcriptomics data
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param pie_scale Scale of each pie to plot. Default is \code{1}.
#' @param xy_ratio Ratio of y and x coordinates. Default is \code{1}.
#' @param color Filled of colors for pie plot, length of \code{color} must be equal to the number of unique cell types in \code{sc_celltype}.
#' @import methods
#' @export

setGeneric("plot_st_pie", def = function(object, pie_scale = 1, xy_ratio = 1, color = NULL) {
    standardGeneric("plot_st_pie")
})

#' @title Plot spatial distribution of gene
#'
#' @description Plot spatial distribution of genes for transcriptomics data
#' @details Please set \code{if_use_newmeta} as \code{FALSE} to plot the spatial distribution of gene before \code{\link{dec_celltype}} for spot-based data.
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param gene Symbol of gene, e.g., 'AKT1'.
#' @param size Point size. Default is \code{1}.
#' @param color_low Color for the lowest value.
#' @param color_mid Color for the middle value for using \code{scale_color_gradient2}. Default is \code{NULL}.
#' @param color_high Color for the highest value.
#' @param color_midpoint Value for the middle scale. Default is \code{NULL}.
#' @param if_use_newmeta Whether to use newmeta to plot the spatial distribution of gene after \code{\link{dec_celltype}} for spot-based data. Default is \code{TRUE}.
#' @param celltype gene in which celltype to plot. Default is \code{NULL}. Set \code{Nif_use_newmeta} TRUE when using this parameter.
#' @param if_plot_others Whether to plot other cells when to use defined \code{celltype}.
#' @import methods
#' @export

setGeneric("plot_st_gene", def = function(object, gene, size = 1, color_low = NULL, color_mid = NULL,
    color_high = NULL, color_midpoint = NULL, if_use_newmeta = T, celltype = NULL, if_plot_others = T) {
    standardGeneric("plot_st_gene")
})

#' @title Plot spatial distribution of a single cell type
#'
#' @description Plot spatial distribution of a single predicted cell types for transcriptomics data
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param celltype Name of cell type in the \code{sc_celltype}.
#' @param size Point size. Default is \code{1}.
#' @param color_celltype Color for the celltype of interest.
#' @param color_others Color for the others.
#' @import methods
#' @export

setGeneric("plot_st_celltype", def = function(object, celltype, size = 1, color_celltype = NULL,
    color_others = NULL) {
    standardGeneric("plot_st_celltype")
})

#' @title Plot spatial density of a single cell type
#'
#' @description Plot spatial density of a single predicted cell types for transcriptomics data
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
#' @import methods
#' @export

setGeneric("plot_st_celltype_density", def = function(object, celltype, type, if_plot_point = T,
    point_color = NULL, point_size = 1, color_low = "grey", color_mid = NULL, color_high = "blue",
    color_midpoint = NULL, size = 1) {
    standardGeneric("plot_st_celltype_density")
})

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
#' @import methods
#' @export

setGeneric("plot_st_celltype_percent", def = function(object, celltype, size = 1, color_low = NULL,
    color_mid = NULL, color_high = NULL, color_midpoint = NULL) {
    standardGeneric("plot_st_celltype_percent")
})

#' @title Plot spatial distribution of all cell types
#'
#' @description Plot spatial distribution of all predicted cell types for transcriptomics data
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param size Point size. Default is \code{1}.
#' @param color Color for all predicted cell types.
#' @import methods
#' @export

setGeneric("plot_st_celltype_all", def = function(object, size = 1, color = NULL) {
    standardGeneric("plot_st_celltype_all")
})

#' @title Plot heatpmap of correlation between marker genes and cell types
#'
#' @description Plot heatpmap of correlation between the expression of marker genes and the predicted score of cell types among all spatial cells or spots.
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param marker_genes A character containing the known marker genes to plot, provide at least two marker genes of interest.
#' @param celltypes A character containing name of cell type in the \code{sc_celltype}. Default is to plot all cell types.
#' @param color_low Color for the lowest value.
#' @param color_mid Color for the middle value for using \code{scale_color_gradient2}. Default is \code{NULL}.
#' @param color_high Color for the highest value.
#' @param if_use_newmeta Whether to use newmeta o plot the spatial distribution of gene after \code{\link{dec_celltype}} for spot-based data. Default is \code{FALSE}.
#' @param scale Character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Corresponding values are 'row', 'column' and 'none'.
#' @param if_show_top Whether to plot a symbol to the highest value across rows or columns. Default is \code{TRUE}.
#' @param top_direction Direction to identify the highest value, select \code{'row'} or \code{'column'}.
#' @param border_color Color of the cell border. Default is \code{'NA'}.
#' @import methods
#' @export

setGeneric("plot_st_cor_heatmap", def = function(object, marker_genes, celltypes, color_low = NULL,
    color_mid = NULL, color_high = NULL, if_use_newmeta = F, scale = "none", if_show_top = T,
    top_direction = "row", border_color = NA) {
    standardGeneric("plot_st_cor_heatmap")
})

#' @title Plot cell-cell distribution
#'
#' @description Plot spatial distribution of celltype_sender and celltype_receiver
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param celltype_sender Name of celltype_sender.
#' @param celltype_receiver Name of celltype_receiver.
#' @param color Color for celltype_sender, celltype_receiver, and others. Three values.
#' @param size Point size. Default is \code{1}.
#' @param if_plot_others Whether to plot others. Default is \code{TRUE}.
#' @param if_plot_density Whether to plot marginal density plots. Default is \code{TRUE}.
#' @param if_plot_edge Whether to plot edge between neighbors. Default is \code{TRUE}.
#' @param arrow_length Arrow length.
#' @param plot_cells Which cells to plot. Default is all cells. Input a character vector of cell names to plot.
#' @import methods
#' @export

setGeneric("plot_ccdist", def = function(object, celltype_sender, celltype_receiver, color = NULL,
    size = 1, if_plot_others = T, if_plot_density = T, if_plot_edge = T, arrow_length = 0.05,
    plot_cells = NULL) {
    standardGeneric("plot_ccdist")
})

#' @title Plot LR pairs
#'
#' @description Plot LR pairs of celltype_sender and celltype_receiver
#' @param object SpaTalk object generated from \code{\link{dec_cci}}.
#' @param celltype_sender Name of celltype_sender.
#' @param celltype_receiver Name of celltype_receiver.
#' @param top_lrpairs Number of top lrpairs for plotting. Default is \code{20}.
#' @param color Color for the cells in heatmap.
#' @param border_color color of cell borders on heatmap, use NA if no border should be drawn.
#' @param type Set 'sig' to plot significant LR pairs or set 'number' to plot the number of spatial LR interactions.
#' @param fontsize_number fontsize of the numbers displayed in cells.
#' @param number_color color of the text.
#' @export

setGeneric("plot_cci_lrpairs", def = function(object, celltype_sender, celltype_receiver, top_lrpairs = 20,
    color = NULL, border_color = "black", type = NULL, fontsize_number = 1, number_color = "black") {
    standardGeneric("plot_cci_lrpairs")
})

#' @title Plot LR pair
#'
#' @description Plot LR pair between celltype_sender and celltype_receiver
#' @param object SpaTalk object generated from \code{\link{dec_cci}}.
#' @param celltype_sender Name of celltype_sender.
#' @param celltype_receiver Name of celltype_receiver.
#' @param ligand Name of ligand from celltype_sender.
#' @param receptor Name of receptor from celltype_receiver.
#' @param color Color for ligand, receptor, and others. Three values.
#' @param size Point size. Default is \code{1}.
#' @param if_plot_density Whether to plot marginal density plots. Default is \code{TRUE}.
#' @param if_plot_edge Whether to plot edge between neighbors. Default is \code{TRUE}.
#' @param arrow_length Arrow length.
#' @param plot_cells Which cells to plot. Default is all cells. Input a character vector of cell names to plot.
#' @import methods
#' @export

setGeneric("plot_lrpair", def = function(object, celltype_sender, celltype_receiver, ligand, receptor,
    color = NULL, size = 1, if_plot_density = T, if_plot_edge = T, arrow_length = 0.05, plot_cells = NULL) {
    standardGeneric("plot_lrpair")
})

#' @title Plot spatial distance of LR pair with vlnplot
#'
#' @description Plot spatial distance of LR pair between expressed senders and receivers and between expressed cell-cell pairs.
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param celltype_sender Name of celltype_sender.
#' @param celltype_receiver Name of celltype_receiver.
#' @param ligand Name of ligand from celltype_sender.
#' @param receptor Name of receptor from celltype_receiver.
#' @param vln_color Color for violins. Two values.
#' @param if_plot_boxplot Whether to plot boxplot. Default is \code{TRUE}.
#' @param box_width Box width. Default is \code{0.2}.
#' @export

setGeneric("plot_lrpair_vln", def = function(object, celltype_sender, celltype_receiver, ligand,
    receptor, vln_color = NULL, if_plot_boxplot = T, box_width = 0.2) {
    standardGeneric("plot_lrpair_vln")
})

#' @title Plot LR and downstream pathways
#'
#' @description Plot LR and downstream pathways
#' @param object SpaTalk object generated from \code{\link{dec_cci}}.
#' @param celltype_sender Name of celltype_sender.
#' @param celltype_receiver Name of celltype_receiver.
#' @param ligand Name of ligand from celltype_sender.
#' @param receptor Name of receptor from celltype_receiver.
#' @param color Color for points Two values.
#' @param size Size of points.
#' @param arrow_length Arrow length.
#' @import methods
#' @export

setGeneric("plot_lr_path", def = function(object, celltype_sender, celltype_receiver, ligand,
    receptor, color = NULL, size = 5, arrow_length = 0.1) {
    standardGeneric("plot_lr_path")
})

#' @title Get LR and downstream pathways
#'
#' @description Get LR and downstream pathways and get p value of receptor-related pathways with LR-target genes by the Fisher-exact test.
#' @param object SpaTalk object generated from \code{\link{dec_cci}}.
#' @param celltype_sender Name of celltype_sender.
#' @param celltype_receiver Name of celltype_receiver.
#' @param ligand Name of ligand from celltype_sender.
#' @param receptor Name of receptor from celltype_receiver.
#' @param min_gene_num Min genes number for each pathway.
#' @return A list containing two data.frame. One is LR and downstream pathways, another is the p value of receptor-related pathways with LR-target genes.
#' @import methods
#' @export

setGeneric("get_lr_path", def = function(object, celltype_sender, celltype_receiver, ligand, receptor,
    min_gene_num = 5) {
    standardGeneric("get_lr_path")
})

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
#' @import methods
#' @export

setGeneric("plot_path2gene", def = function(object, celltype_sender, celltype_receiver,
  ligand, receptor, min_gene_num = 5, pvalue = 0.5, color = NULL, color_flow = "blue") {
  standardGeneric("plot_path2gene")
})
