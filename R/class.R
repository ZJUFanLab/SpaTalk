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
    cellpair = "list", dist = "matrix", lrpair = "data.frame", tf = "data.frame",
    lr_path = "list"), prototype(data = list(), meta = list(), para = list(), coef = matrix(),
    cellpair = list(), dist = matrix(), lrpair = data.frame(), tf = data.frame(),
    lr_path = list()))
