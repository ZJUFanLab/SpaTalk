#' @title Generate pseudo spot st_data
#'
#' @description Generate pseudo spot st_data with single-cell st_data
#' @param st_data A data.frame or matrix or dgCMatrix containing counts of spatial transcriptomics, each column representing a cell, each row representing a gene.
#' @param st_meta A data.frame containing coordinate of spatial transcriptomics with three columns, \code{'cell'}, \code{'x'}, \code{'y'}, and \code{celltype}.
#' @param x_min Min value of x axis.
#' @param x_res Resolution of x coordinate.
#' @param x_max Max value of x axis.
#' @param y_min Min value of y axis.
#' @param y_res Resolution of y coordinate.
#' @param y_max Max value of y axis.
#' @return A list of spot st_data and st_meta
#' @export
#' @import Matrix
#' @importFrom reshape2 dcast
#' @importFrom methods as

generate_spot <- function(st_data, st_meta, x_min, x_res, x_max, y_min, y_res, y_max) {
    if (is(st_data, "data.frame")) {
        st_data <- methods::as(as.matrix(st_data), "dgCMatrix")
    }
    if (is(st_data, "matrix")) {
        st_data <- methods::as(st_data, "dgCMatrix")
    }
    if (!is(st_data, "dgCMatrix")) {
        stop("st_data must be a data.frame or matrix or dgCMatrix!")
    }
    if (is(st_meta, "data.frame")) {
        if (!all(c("cell","x","y") %in% colnames(st_meta))) {
            stop("Please provide a correct st_meta data.frame! See demo_st_sc_meta()!")
        }
    } else {
        stop("st_meta must be a data.frame!")
    }
    x_range <- seq(from = x_min, to = x_max, by = x_res)
    y_range <- seq(from = y_min, to = y_max, by = y_res)
    # x and y resolution to construct spot data
    test_spot_plot <- data.frame(spot = "x", x = 0, y = 0, stringsAsFactors = F)
    for (i in 1:(length(x_range) - 1)) {
        x1 <- paste0("x", rep(i, length(y_range) - 1))
        test_spot_plot1 <- data.frame(spot = x1, x = x_range[i] + x_res, y = 0, stringsAsFactors = F)
        for (j in 1:(length(y_range) - 1)) {
            y1 <- paste0(test_spot_plot1$spot[j], "_", "y", j)
            test_spot_plot1$spot[j] <- y1
            test_spot_plot1$y[j] <- y_range[j] + y_res
        }
        test_spot_plot <- rbind(test_spot_plot, test_spot_plot1)
    }
    test_spot_plot <- test_spot_plot[-1, ]
    for (i in 1:nrow(st_meta)) {
        test_data_x <- max(which(x_range <= st_meta$x[i]))
        test_data_y <- max(which(y_range <= st_meta$y[i]))
        st_meta$spot[i] <- paste0("x", test_data_x, "_", "y", test_data_y)
    }
    test_spot_meta <- as.data.frame(table(st_meta$spot), stringsAsFactors = F)
    test_spot_meta1 <- reshape2::dcast(data = st_meta[, c("spot", "celltype")], formula = spot ~ celltype)
    test_spot_meta <- cbind(test_spot_meta, test_spot_meta1[, -1])
    # test_spot_data -- sum
    test_spot_data <- list()
    for (i in 1:nrow(test_spot_meta)) {
        test_spot_cell <- st_data[, st_meta[st_meta$spot == test_spot_meta$Var1[i], ]$cell]
        if (is(test_spot_cell, "dgCMatrix")) {
            test_spot_sum <- rowSums(test_spot_cell)
        } else {
            test_spot_sum <- test_spot_cell
        }
        test_spot_data[[i]] <- test_spot_sum
        names(test_spot_data)[i] <- test_spot_meta$Var1[i]
    }
    test_spot_data <- as.data.frame(test_spot_data, stringsAsFactors = F)
    test_spot_data <- as(as.matrix(test_spot_data), "dgCMatrix")
    # generate x and y
    rownames(test_spot_plot) <- test_spot_plot$spot
    test_spot_plot <- test_spot_plot[test_spot_meta$Var1, ]
    test_spot_real <- test_spot_meta
    colnames(test_spot_real)[c(1, 2)] <- c("spot", "cell_real")
    test_spot_meta <- test_spot_plot
    test_spot_meta <- cbind(test_spot_meta, test_spot_real[, -1])
    return(list(st_data = test_spot_data, st_meta = test_spot_meta))
}

#' @title SpaTalk object
#'
#' @description create SpaTalk object using spatial transcriptomics data.
#' @param st_data A data.frame or matrix or dgCMatrix containing counts of spatial transcriptomics, each column representing a spot or a cell, each row representing a gene.
#' @param st_meta A data.frame containing coordinate of spatial transcriptomics with three columns, namely \code{'spot'}, \code{'x'}, \code{'y'} for spot-based spatial transcriptomics data or \code{'cell'}, \code{'x'}, \code{'y'} for single-cell spatial transcriptomics data.
#' @param species A character meaning species of the spatial transcriptomics data.\code{'Human'} or \code{'Mouse'}.
#' @param if_st_is_sc A logical meaning if it is single-cell spatial transcriptomics data. \code{TRUE} is \code{FALSE}.
#' @param spot_max_cell A integer meaning max cell number for each plot to predict. If \code{if_st_sc} is \code{FALSE}, please determine the \code{spot_max_cell}. For 10X (55um), we recommend 30. For Slide-seq, we recommend 1.
#' @param celltype A character containing the cell type of ST data. To skip the deconvolution step and directly infer cell-cell communication, please define the cell type. Default is `NULL`.
#' @return SpaTalk object
#' @importFrom  methods as
#' @import Matrix
#' @export

createSpaTalk <- function(st_data, st_meta, species, if_st_is_sc, spot_max_cell, celltype = NULL) {
    if (is(st_data, "data.frame")) {
        st_data <- methods::as(as.matrix(st_data), "dgCMatrix")
    }
    if (is(st_data, "matrix")) {
        st_data <- methods::as(st_data, "dgCMatrix")
    }
    if (!is(st_data, "dgCMatrix")) {
        stop("st_data must be a data.frame or matrix or dgCMatrix!")
    }
    if (!is.data.frame(st_meta)) {
        stop("st_meta is not a data frame!")
    }
    # check st_data and st_meta
    if (if_st_is_sc) {
        if (!all(c("cell", "x", "y") == colnames(st_meta))) {
            stop("Please provide a correct st_meta data.frame! See demo_st_sc_meta()!")
        }
        if (!all(colnames(st_data) == st_meta$cell)) {
            stop("colnames(st_data) is not consistent with st_meta$cell!")
        }
        st_type <- "single-cell"
        spot_max_cell <- 1
    } else {
        if (!all(c("spot", "x", "y") == colnames(st_meta))) {
            stop("Please provide a correct st_meta data.frame! See demo_st_meta()!")
        }
        if (!all(colnames(st_data) == st_meta$spot)) {
            stop("colnames(st_data) is not consistent with st_meta$spot!")
        }
        st_type <- "spot"
    }
    if (is.null(spot_max_cell)) {
        stop("Please provide the spot_max_cell!")
    }
    st_data <- st_data[which(rowSums(st_data) > 0), ]
    if (nrow(st_data) == 0) {
        stop("No expressed genes in st_data!")
    }
    # st_meta
    st_meta$nFeatures <- as.numeric(apply(st_data, 2, .percent_cell))
    st_meta$label <- "-"
    st_meta$cell_num <- spot_max_cell
    st_meta$celltype <- "unsure"
    if (!is.null(celltype)) {
        if (!is.character(celltype)) {
            stop("celltype must be a character with length equal to ST data!")
        }
        if (length(celltype) != nrow(st_meta)) {
            stop("Length of celltype must be equal to nrow(st_meta)!")
        }
        celltype_new <- .rename_chr(celltype)
        warning_info <- .show_warning(celltype, celltype_new)
        if (!is.null(warning_info)) {
            warning(warning_info)
        }
        st_meta$celltype <- celltype_new
        if_skip_dec_celltype <- TRUE
    } else {
        if_skip_dec_celltype <- FALSE
    }
    st_meta[, 1] <- .rename_chr(st_meta[, 1])
    colnames(st_data) <- st_meta[,1]
    # generate SpaTalk object
    object <- new("SpaTalk", data = list(rawdata = st_data), meta = list(rawmeta = st_meta),
        para = list(species = species, st_type = st_type, spot_max_cell = spot_max_cell, if_skip_dec_celltype = if_skip_dec_celltype))
    return(object)
}

#' @title Decomposing cell type for spatial transcriptomics data
#'
#' @description Identify the cellular composition for single-cell or spot-based spatial transcriptomics data with non-negative regression.
#' @param object SpaTalk object generated from \code{\link{createSpaTalk}}.
#' @param sc_data A A data.frame or matrix or dgCMatrix containing counts of single-cell RNA-seq data as the reference, each column representing a cell, each row representing a gene.
#' @param sc_celltype A character containing the cell type of the reference single-cell RNA-seq data.
#' @param min_percent Min percent to predict new cell type for single-cell st_data or predict new cell for spot-based st_data. Default is \code{0.5}.
#' @param min_nFeatures Min number of expressed features/genes for each spot/cell in \code{st_data}. Default is \code{10}.
#' @param if_use_normalize_data Whether to use normalized \code{st_data} and \code{sc_data} with Seurat normalization. Default is \code{TRUE}. set it \code{FALSE} when the st_data and sc_data are already normalized matrix with other methods.
#' @param if_use_hvg Whether to use highly variable genes for non-negative regression. Default is \code{FALSE}.
#' @param if_retain_other_genes Whether to retain other genes which are not overlapped between sc_data and st_data when reconstructing the single-cell ST data. Default is \code{FALSE}. Set it \code{TRUE} to obtain the constructed single-cell ST data with genes consistent with that in sc_data.
#' @param if_doParallel Use doParallel. Default is TRUE.
#' @param use_n_cores Number of CPU cores to use. Default is all cores - 2.
#' @param iter_num Number of iteration to generate the single-cell data for spot-based data. Default is \code{1000}.
#' @param method 1 means using the SpaTalk deconvolution method, 2 means using RCTD, 3 means using Seurat, 4 means using SPOTlight, 5 means using deconvSeq, 6 means using stereoscope, 7 means using cell2location
#' @param env When method set to 6, namely use stereoscope python package to deconvolute, please define the python environment of installed stereoscope. Default is the 'base' environment. Anaconda is recommended. When method set to 7, namely use cell2location python package to deconvolute, please install cell2location to "base" environment.
#' @param anaconda_path When using stereoscope, please define the \code{env} parameter as well as the path to anaconda. Default is "~/anaconda3"
#' @param dec_result A matrix of deconvolution result from other upcoming methods, row represents spots or cells, column represents cell types of scRNA-seq reference. See \code{\link{demo_dec_result}}
#' @return SpaTalk object containing the decomposing results.
#' @import Matrix progress methods Seurat foreach doParallel parallel iterators readr
#' @importFrom crayon cyan green
#' @importFrom stringr str_replace_all
#' @importFrom NNLM nnlm
#' @importFrom stats dist
#' @export

dec_celltype <- function(object, sc_data, sc_celltype, min_percent = 0.5, min_nFeatures = 10, if_use_normalize_data = T, if_use_hvg = F,
    if_retain_other_genes = F, if_doParallel = T, use_n_cores = NULL, iter_num = 1000, method = 1, env = "base", anaconda_path = "~/anaconda3", dec_result = NULL) {
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    if_skip_dec_celltype <- object@para$if_skip_dec_celltype
    if (if_skip_dec_celltype) {
        stop("Do not perform dec_celltype() when providing celltype in createSpaTalk()!")
    }
    if (if_doParallel) {
        if (is.null(use_n_cores)) {
            n_cores <- parallel::detectCores()
            n_cores <- floor(n_cores/4)
        } else {
            n_cores <- use_n_cores
        }
        n.threads <- n_cores
        if (n_cores < 2) {
            if_doParallel <- F
            n.threads <- 0
        }
    } else {
        n_cores <- 1
        n.threads <- 0
    }
    st_data <- object@data[["rawdata"]]
    st_meta <- object@meta[["rawmeta"]]
    # dist
    st_dist <- .st_dist(st_meta)
    if (min(st_meta$nFeatures) < min_nFeatures) {
        st_meta[st_meta$nFeatures < min_nFeatures, ]$label <- "less nFeatures"
    }
    st_type <- object@para[["st_type"]]
    if (is(sc_data, "data.frame")) {
        sc_data <- methods::as(as.matrix(sc_data), "dgCMatrix")
    }
    if (is(sc_data, "matrix")) {
        sc_data <- methods::as(sc_data, "dgCMatrix")
    }
    if (!is(sc_data, "dgCMatrix")) {
        stop("sc_data must be a data.frame or matrix or dgCMatrix!")
    }
    if (!is.character(sc_celltype)) {
        stop("sc_celltype is not a character!")
    }
    if (ncol(sc_data) != length(sc_celltype)) {
        stop("ncol(sc_data) is not consistent with length(sc_celltype)!")
    }
    if (min_percent >= 1 | min_percent <= 0) {
        stop("Please provide a correct min_percent, ranging from 0 to 1!")
    }
    sc_data <- sc_data[which(rowSums(sc_data) > 0), ]
    if (nrow(sc_data) == 0) {
        stop("No expressed genes in sc_data!")
    }
    colnames(sc_data) <- .rename_chr(colnames(sc_data))
    sc_celltype_new <- .rename_chr(sc_celltype)
    warning_info <- .show_warning(sc_celltype, sc_celltype_new)
    if (!is.null(warning_info)) {
        warning(warning_info)
    }
    sc_celltype <- data.frame(cell = colnames(sc_data), celltype = sc_celltype_new, stringsAsFactors = F)
    object@data$rawndata <- st_data
    if (if_use_normalize_data == T) {
        st_ndata <- .normalize_data(st_data)
        sc_ndata <- .normalize_data(sc_data)
    } else {
        st_ndata <- st_data
        sc_ndata <- sc_data
    }
    sc_ndata_raw <- sc_ndata
    genename <- intersect(rownames(st_data), rownames(sc_data))
    if (length(genename) == 0) {
        stop("No overlapped genes between st_data and sc_data!")
    }
    if (length(genename) == 1) {
        stop("Only 1 gene overlapped between st_data and sc_data!")
    }
    if (method == 1 & is.null(dec_result)) {
        object@data$rawndata <- st_ndata
        st_ndata <- st_ndata[genename, ]
        sc_ndata <- sc_ndata[genename, ]
        # performing Non-negative regression
        if (st_type == "single-cell") {
            cat(crayon::cyan("Performing Non-negative regression for each cell", "\n"))
        } else {
            cat(crayon::cyan("Performing Non-negative regression for each spot", "\n"))
        }
        # use hvg
        if (if_use_hvg) {
            sc_hvg <- .hvg(sc_ndata, sc_celltype)
            if (length(sc_hvg) > 0) {
                st_ndata <- st_ndata[sc_hvg, ]
                sc_ndata <- sc_ndata[sc_hvg, ]
            }
        }
        if (if_doParallel) {
            cl <- parallel::makeCluster(n_cores)
            doParallel::registerDoParallel(cl)
        }
        # generate sc_ref data
        sc_ref <- .generate_ref(sc_ndata, sc_celltype, if_doParallel)
        if (if_doParallel) {
            doParallel::stopImplicitCluster()
            parallel::stopCluster(cl)
        }
        # deconvolution
        st_nnlm <- NNLM::nnlm(x = as.matrix(sc_ref), y = as.matrix(st_ndata), method = "lee", loss = "mkl", n.threads = n.threads)
        st_coef <- t(st_nnlm$coefficients)
    }
    if (method == 2 & is.null(dec_result)) {
        cat(crayon::cyan("Using RCTD to deconvolute", "\n"))
        require(spacexr)
        st_coef <- .run_rctd(st_data, sc_data, st_meta, sc_celltype)
    }
    if (method == 3 & is.null(dec_result)) {
        cat(crayon::cyan("Using Seurat to deconvolute", "\n"))
        require(Seurat)
        st_coef <- .run_seurat(st_data, sc_data, sc_celltype)
    }
    if (method == 4 & is.null(dec_result)) {
        cat(crayon::cyan("Using SPOTlight to deconvolute", "\n"))
        require(SPOTlight)
        st_coef <- .run_spotlight(st_data, sc_data, sc_celltype)
    }
    if (method == 5 & is.null(dec_result)) {
        cat(crayon::cyan("Using deconvSeq to deconvolute", "\n"))
        require(deconvSeq)
        st_coef <- .run_deconvSeq(st_data, sc_data, sc_celltype)
    }
    if (method == 6 & is.null(dec_result)) {
        cat(crayon::cyan("Using stereoscope to deconvolute, please install the stereoscope (python package) first!", "\n"))
        # python
        st_coef <- .run_stereoscope(st_data, sc_data, sc_celltype, env, anaconda_path)
    }
    if (method == 7 & is.null(dec_result)) {
        cat(crayon::cyan("Using cell2location to deconvolute, please install the cell2location (python package) first!", "\n"))
        # python
        require(anndata)
        require(Seurat)
        require(reticulate)
        require(sceasy)
        st_coef <- .run_cell2location(st_data, st_meta, sc_data, sc_celltype, env)
    }
    if (!is.null(dec_result)) {
        if (!is.matrix(dec_result)) {
            stop("Please provide a correct dec_result matrix! See demo_dec_result()!")
        }
        dec_colname <- colnames(dec_result)
        dec_colname <- .rename_chr(dec_colname)
        colnames(dec_result) <- dec_colname
        if (!all(dec_colname %in% unique(sc_celltype$celltype))) {
            stop("Celltype name in dec_result must be consistent with the names in scRNA-seq reference!")
        }
        dec_rowname <- rownames(dec_result)
        dec_rowname <- .rename_chr(dec_rowname)
        rownames(dec_result) <- dec_rowname
        if (!all(colnames(st_data) == dec_rowname)) {
            stop("Spot/cell name in dec_result must be consistent with the names in st_meta!")
        }
        st_coef <- dec_result
    }
    coef_name <- colnames(st_coef)
    coef_name <- coef_name[order(coef_name)]
    st_coef <- st_coef[ ,coef_name]
    object@coef <- st_coef
    st_meta <- cbind(st_meta, .coef_nor(st_coef))
    st_ndata <- object@data$rawndata
    sc_ndata <- sc_ndata_raw
    st_ndata <- st_ndata[genename, ]
    sc_ndata <- sc_ndata[genename, ]
    if (st_type == "single-cell") {
        object@meta$rawmeta <- .determine_celltype(st_meta, min_percent)
    } else {
        if (if_doParallel) {
            cl <- parallel::makeCluster(n_cores)
            doParallel::registerDoParallel(cl)
            newmeta <- .generate_newmeta_doParallel(st_meta, st_dist, min_percent)
            doParallel::stopImplicitCluster()
            parallel::stopCluster(cl)
        } else {
            newmeta <- .generate_newmeta(st_meta, st_dist, min_percent)
        }
        if (nrow(newmeta) > 0) {
            newmeta_cell <- .generate_newmeta_cell(newmeta, st_ndata, sc_ndata, sc_celltype, iter_num, n_cores, if_doParallel)
            if (if_retain_other_genes) {
                sc_ndata <- sc_ndata_raw
            } else {
                if (if_use_hvg) {
                    sc_ndata <- sc_ndata_raw
                    sc_ndata <- sc_ndata[genename, ]
                }
            }
            newdata <- sc_ndata[, newmeta_cell$cell_id]
            colnames(newdata) <- newmeta_cell$cell
            object@data$newdata <- methods::as(newdata, Class = "dgCMatrix")
            object@meta$newmeta <- newmeta_cell
            st_meta[st_meta$spot %in% newmeta_cell$spot, ]$celltype <- "sure"
            st_dist <- .st_dist(newmeta_cell)
            object@meta$rawmeta <- st_meta
        } else {
            warning("No new data are generated for the min_percent!")
        }
    }
    cat(crayon::green("***Done***", "\n"))
    object@dist <- st_dist
    object@para$min_percent <- min_percent
    object@para$min_nFeatures <- min_nFeatures
    object@para$if_use_normalize_data <- if_use_normalize_data
    object@para$if_use_hvg <- if_use_hvg
    object@para$iter_num <- iter_num
    return(object)
}

#' @title  Set the expected cell
#'
#' @description Set the expected cell in SpaTalk object
#' @param object SpaTalk object
#' @param value Th number of expected cell for each spot, must be equal to the spot number.
#' @export set_expected_cell
#' @return SpaTalk object

set_expected_cell <- function(object, value) {
  if (!is(object, "SpaTalk")) {
      stop("Invalid class for object: must be 'SpaTalk'!")
  }
  st_meta <- object@meta$rawmeta
  if (length(value) != nrow(st_meta)) {
      stop("The value must be equal to the spot number!")
  }
  object@meta$rawmeta$cell_num <- value
  return(object)
}

#' @title Find lrpairs and pathways
#'
#' @description Find \code{lrpairs} and \code{pathways} with receptors having downstream targets and transcriptional factors.
#' @param object SpaTalk object generated from \code{\link{dec_celltype}}.
#' @param lrpairs A data.frame of the system data containing ligand-receptor pairs of \code{'Human'} and \code{'Mouse'} from CellTalkDB.
#' @param pathways A data.frame of the system data containing gene-gene interactions and pathways from KEGG and Reactome as well as the information of transcriptional factors.
#' @param max_hop Max hop from the receptor to the downstream target transcriptional factor to find for receiving cells. Default is \code{3} for human and \code{4} for mouse.
#' @param if_doParallel Use doParallel. Default is TRUE.
#' @param use_n_cores Number of CPU cores to use. Default is all cores - 2.
#' @return SpaTalk object containing the filtered lrpairs and pathways.
#' @import Matrix progress methods
#' @importFrom crayon cyan green
#' @export find_lr_path

find_lr_path <- function(object, lrpairs, pathways, max_hop = NULL, if_doParallel = T, use_n_cores = NULL) {
    # check input data
    cat(crayon::cyan("Checking input data", "\n"))
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    if (!all(c("ligand", "receptor", "species") %in% names(lrpairs))) {
        stop("Please provide a correct lrpairs data.frame! See demo_lrpairs()!")
    }
    if (!all(c("src", "dest", "pathway", "source", "type", "src_tf", "dest_tf", "species") %in%
        names(pathways))) {
        stop("Please provide a correct pathways data.frame! See demo_pathways()!")
    }
    if (is.null(use_n_cores)) {
        n_cores <- parallel::detectCores()
        n_cores <- floor(n_cores/4)
        if (n_cores < 2) {
            if_doParallel <- FALSE
        }
    } else {
        n_cores <- use_n_cores
    }
    if (if_doParallel) {
        cl <- parallel::makeCluster(n_cores)
        doParallel::registerDoParallel(cl)
    }
    st_data <- .get_st_data(object)
    species <- object@para$species
    lrpair <- lrpairs[lrpairs$species == species, ]
    lrpair <- lrpair[lrpair$ligand %in% rownames(st_data) & lrpair$receptor %in% rownames(st_data), ]
    if (nrow(lrpair) == 0) {
        stop("No ligand-recepotor pairs found in st_data!")
    }
    pathways <- pathways[pathways$species == species, ]
    pathways <- pathways[pathways$src %in% rownames(st_data) & pathways$dest %in% rownames(st_data), ]
    ggi_tf <- pathways[, c("src", "dest", "src_tf", "dest_tf")]
    ggi_tf <- unique(ggi_tf)
    lrpair <- lrpair[lrpair$receptor %in% ggi_tf$src, ]
    if (nrow(lrpair) == 0) {
        stop("No downstream target genes found for receptors!")
    }
    cat(crayon::cyan("Begin to filter lrpairs and pathways", "\n"))
    if (is.null(max_hop)) {
        if (species == "Mouse") {
            max_hop <- 4
        } else {
            max_hop <- 3
        }
    }
    ### find receptor-tf
    res_ggi <- NULL
    receptor_name <- unique(lrpair$receptor)
    if (if_doParallel) {
        res_ggi <- foreach::foreach(i=1:length(receptor_name), .combine = "c", .packages = "Matrix") %dopar% {
            ggi_res <- NULL
            lr_receptor <- receptor_name[i]
            ggi_tf1 <- ggi_tf[ggi_tf$src == lr_receptor, ]
            ggi_tf1 <- unique(ggi_tf1[ggi_tf1$dest %in% rownames(st_data), ])
            if (nrow(ggi_tf1) > 0) {
                n <- 0
                ggi_tf1$hop <- n + 1
                while (n <= max_hop) {
                    ggi_res <- rbind(ggi_res, ggi_tf1)
                    ggi_tf1 <- ggi_tf[ggi_tf$src %in% ggi_tf1$dest, ]
                    ggi_tf1 <- unique(ggi_tf1[ggi_tf1$dest %in% rownames(st_data), ])
                    if (nrow(ggi_tf1) == 0) {
                        break
                    }
                    ggi_tf1$hop <- n + 2
                    n <- n + 1
                }
                ggi_res <- unique(ggi_res)
                ggi_res_yes <- ggi_res[ggi_res$src_tf == "YES" | ggi_res$dest_tf == "YES", ]
                if (nrow(ggi_res_yes) > 0) {
                    lr_receptor
                }
            }
        }
    } else {
        for (i in 1:length(receptor_name)) {
            ggi_res <- NULL
            lr_receptor <- receptor_name[i]
            ggi_tf1 <- ggi_tf[ggi_tf$src == lr_receptor, ]
            ggi_tf1 <- unique(ggi_tf1[ggi_tf1$dest %in% rownames(st_data), ])
            if (nrow(ggi_tf1) > 0) {
                n <- 0
                ggi_tf1$hop <- n + 1
                while (n <= max_hop) {
                    ggi_res <- rbind(ggi_res, ggi_tf1)
                    ggi_tf1 <- ggi_tf[ggi_tf$src %in% ggi_tf1$dest, ]
                    ggi_tf1 <- unique(ggi_tf1[ggi_tf1$dest %in% rownames(st_data), ])
                    if (nrow(ggi_tf1) == 0) {
                        break
                    }
                    ggi_tf1$hop <- n + 2
                    n <- n + 1
                }
                ggi_res <- unique(ggi_res)
                ggi_res_yes <- ggi_res[ggi_res$src_tf == "YES" | ggi_res$dest_tf == "YES", ]
                if (nrow(ggi_res_yes) > 0) {
                    res_ggi <- c(res_ggi, lr_receptor)
                }
            }
        }
    }
    if (length(res_ggi) == 0) {
        stop("No downstream transcriptional factors found for receptors!")
    }
    cat(crayon::green("***Done***", "\n"))
    if (if_doParallel) {
        doParallel::stopImplicitCluster()
        parallel::stopCluster(cl)
    }
    lrpair <- lrpair[lrpair$receptor %in% res_ggi, ]
    if (nrow(lrpair) == 0) {
        stop("No ligand-recepotor pairs found!")
    }
    object@lr_path <- list(lrpairs = lrpair, pathways = pathways)
    object@para$max_hop <- max_hop
    if_skip_dec_celltype <- object@para$if_skip_dec_celltype
    if (if_skip_dec_celltype) {
        st_meta <- object@meta$rawmeta
        st_dist <- .st_dist(st_meta)
        object@dist <- st_dist
    }
    return(object)
}

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
#' @param if_doParallel Use doParallel. Default is TRUE.
#' @param use_n_cores Number of CPU cores to use. Default is all cores - 2.
#' @return SpaTalk object containing the inferred LR pairs and pathways.
#' @import Matrix progress methods
#' @importFrom crayon cyan green
#' @export

dec_cci <- function(object, celltype_sender, celltype_receiver, n_neighbor = 10, min_pairs = 5,
    min_pairs_ratio = 0, per_num = 1000, pvalue = 0.05, co_exp_ratio = 0.1, if_doParallel = T, use_n_cores = NULL) {
    # check input data
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    if (is.null(use_n_cores)) {
        n_cores <- parallel::detectCores()
        n_cores <- floor(n_cores/4)
        if (n_cores < 2) {
            if_doParallel <- FALSE
        }
    } else {
        n_cores <- use_n_cores
    }
    if (if_doParallel) {
        cl <- parallel::makeCluster(n_cores)
        doParallel::registerDoParallel(cl)
    }
    st_meta <- .get_st_meta(object)
    st_data <- .get_st_data(object)
    celltype_dist <- object@dist
    if (!celltype_sender %in% st_meta$celltype) {
        stop("Please provide a correct celltype_sender")
    }
    if (!celltype_receiver %in% st_meta$celltype) {
        stop("Please provide a correct celltype_receiver")
    }
    cell_pair <- .get_cellpair(celltype_dist, st_meta, celltype_sender, celltype_receiver, n_neighbor)
    cell_sender <- st_meta[st_meta$celltype == celltype_sender, ]
    cell_receiver <- st_meta[st_meta$celltype == celltype_receiver, ]
    cell_pair_all <- nrow(cell_sender) * nrow(cell_receiver)/2
    if (nrow(cell_pair) <= min_pairs) {
        stop(paste0("Cell pairs found between ", celltype_sender, " and ", celltype_receiver, " less than min_pairs!"))
    }
    if (nrow(cell_pair) <= cell_pair_all * min_pairs_ratio) {
        stop(paste0("Cell pairs found between ", celltype_sender, " and ", celltype_receiver, " less than min_pairs_ratio!"))
    }
    lrdb <- object@lr_path$lrpairs
    pathways <- object@lr_path$pathways
    ggi_tf <- unique(pathways[, c("src", "dest", "src_tf", "dest_tf")])
    cat(crayon::cyan("Begin to find LR pairs", "\n"))
    ### [1] LR distance
    if (if_doParallel) {
        lrdb <- .lr_distance_doParallel(st_data, cell_pair, lrdb, celltype_sender, celltype_receiver, per_num, pvalue)
    } else {
        lrdb <- .lr_distance(st_data, cell_pair, lrdb, celltype_sender, celltype_receiver, per_num, pvalue)
    }
    ### [2] Downstream targets and TFs
    max_hop <- object@para$max_hop
    receptor_tf <- NULL
    if (nrow(lrdb) > 0) {
        if (if_doParallel) {
            receptor_tf <- .get_tf_res_doParallel(celltype_sender, celltype_receiver, lrdb, ggi_tf, cell_pair, st_data, max_hop, co_exp_ratio)
        } else {
            receptor_tf <- .get_tf_res(celltype_sender, celltype_receiver, lrdb, ggi_tf, cell_pair, st_data, max_hop, co_exp_ratio)
        }
        if (is.null(receptor_tf)) {
            stop(paste0("No LR pairs found between ", celltype_sender, " and ", celltype_receiver))
        }
        # calculate score
        lrdb <- .get_score(lrdb, receptor_tf)
    } else {
        stop(paste0("No LR pairs found between ", celltype_sender, " and ", celltype_receiver))
    }
    lrpair <- object@lrpair
    if (nrow(lrpair) == 0) {
        object@lrpair <- lrdb
    } else {
        lrpair <- rbind(lrpair, lrdb)
        object@lrpair <- unique(lrpair)
    }
    tf <- object@tf
    if (nrow(tf) == 0) {
        object@tf <- receptor_tf
    } else {
        tf <- rbind(tf, receptor_tf)
        object@tf <- unique(tf)
    }
    object@cellpair[[paste0(celltype_sender, " -- ", celltype_receiver)]] <- cell_pair
    if (if_doParallel) {
        doParallel::stopImplicitCluster()
        parallel::stopCluster(cl)
    }
    object@para$min_pairs <- min_pairs
    object@para$min_pairs_ratio <- min_pairs_ratio
    object@para$per_num <- per_num
    object@para$pvalue <- pvalue
    object@para$co_exp_ratio <- co_exp_ratio
    return(object)
}

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
#' @param if_doParallel Use doParallel. Default is TRUE.
#' @param use_n_cores Number of CPU cores to use. Default is all cores - 2.
#' @return SpaTalk object containing the inferred LR pairs and pathways.
#' @import Matrix progress methods
#' @importFrom crayon cyan combine_styles
#' @export

dec_cci_all <- function(object, n_neighbor = 10, min_pairs = 5, min_pairs_ratio = 0,
    per_num = 1000, pvalue = 0.05, co_exp_ratio = 0.1, if_doParallel = T, use_n_cores = NULL) {
    # check input data
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    if (is.null(use_n_cores)) {
        n_cores <- parallel::detectCores()
        n_cores <- floor(n_cores/4)
        if (n_cores < 2) {
            if_doParallel <- FALSE
        }
    } else {
        n_cores <- use_n_cores
    }
    if (if_doParallel) {
        cl <- parallel::makeCluster(n_cores)
        doParallel::registerDoParallel(cl)
    }
    st_meta <- .get_st_meta(object)
    st_data <- .get_st_data(object)
    celltype_dist <- object@dist
    # generate pair-wise cell types
    cellname <- unique(st_meta$celltype)
    celltype_pair <- NULL
    for (i in 1:length(cellname)) {
        d1 <- data.frame(celltype_sender = rep(cellname[i], length(cellname)), celltype_receiver = cellname,
            stringsAsFactors = F)
        celltype_pair <- rbind(celltype_pair, d1)
    }
    celltype_pair <- celltype_pair[celltype_pair$celltype_sender != celltype_pair$celltype_receiver, ]
    cat(paste0("Note: there are ", length(cellname), " cell types and ", nrow(celltype_pair), " pair-wise cell pairs"), "\n")
    pathways <- object@lr_path$pathways
    ggi_tf <- unique(pathways[, c("src", "dest", "src_tf", "dest_tf")])
    cat(crayon::cyan("Begin to find LR pairs", "\n"))
    if (if_doParallel) {
        all_res <- foreach::foreach(i=1:nrow(celltype_pair), .packages = c("Matrix", "reshape2"),
            .export = c(".get_cellpair", ".lr_distance", ".get_tf_res", ".get_score")) %dopar% {
            celltype_sender <- celltype_pair$celltype_sender[i]
            celltype_receiver <- celltype_pair$celltype_receiver[i]
            cell_pair <- .get_cellpair(celltype_dist, st_meta, celltype_sender, celltype_receiver, n_neighbor)
            cell_sender <- st_meta[st_meta$celltype == celltype_sender, ]
            cell_receiver <- st_meta[st_meta$celltype == celltype_receiver, ]
            cell_pair_all <- nrow(cell_sender) * nrow(cell_receiver)/2
            if (nrow(cell_pair) > min_pairs) {
                if (nrow(cell_pair) > cell_pair_all * min_pairs_ratio) {
                    lrdb <- object@lr_path$lrpairs
                    ### [1] LR distance
                    lrdb <- .lr_distance(st_data, cell_pair, lrdb, celltype_sender, celltype_receiver, per_num, pvalue)
                    ### [2] Downstream targets and TFs
                    max_hop <- object@para$max_hop
                    receptor_tf <- NULL
                    if (nrow(lrdb) > 0) {
                        receptor_tf <- .get_tf_res(celltype_sender, celltype_receiver, lrdb, ggi_tf, cell_pair, st_data, max_hop, co_exp_ratio)
                        if (!is.null(receptor_tf)) {
                            # calculate score
                            lrdb <- .get_score(lrdb, receptor_tf)
                        } else {
                            lrdb <- NULL
                        }
                    }
                    if (is.data.frame(lrdb)) {
                        if (nrow(lrdb) > 0) {
                            list(lrdb = lrdb, receptor_tf=receptor_tf, cell_pair=cell_pair)
                        }
                    }
                }
            }
        }
        doParallel::stopImplicitCluster()
        parallel::stopCluster(cl)
        res_receptor_tf <- NULL
        res_lrpair <- NULL
        res_cellpair <- list()
        m <- 0
        for (i in 1:length(all_res)) {
            all_res1 <- all_res[[i]]
            if (!is.null(all_res1)) {
                m <- m+1
                res_lrpair <- rbind(res_lrpair, all_res1[[1]])
                res_receptor_tf <- rbind(res_receptor_tf, all_res1[[2]])
                res_cellpair[[m]] <- all_res1[[3]]
                names(res_cellpair)[m] <- paste0(unique(all_res1[[1]]$celltype_sender), " -- ", unique(all_res1[[1]]$celltype_receiver))
            }
        }
        if (!is.null(res_lrpair)) {
            object@lrpair <- res_lrpair
        }
        if (!is.null(res_receptor_tf)) {
            object@tf <- res_receptor_tf
        }
        if (length(res_cellpair) > 0) {
            object@cellpair <- res_cellpair
        }
    } else {
        for (i in 1:nrow(celltype_pair)) {
            celltype_sender <- celltype_pair$celltype_sender[i]
            celltype_receiver <- celltype_pair$celltype_receiver[i]
            cell_pair <- .get_cellpair(celltype_dist, st_meta, celltype_sender, celltype_receiver, n_neighbor)
            cell_sender <- st_meta[st_meta$celltype == celltype_sender, ]
            cell_receiver <- st_meta[st_meta$celltype == celltype_receiver, ]
            cell_pair_all <- nrow(cell_sender) * nrow(cell_receiver)/2
            if (nrow(cell_pair) > min_pairs) {
                if (nrow(cell_pair) > cell_pair_all * min_pairs_ratio) {
                    lrdb <- object@lr_path$lrpairs
                    ### [1] LR distance
                    lrdb <- .lr_distance(st_data, cell_pair, lrdb, celltype_sender, celltype_receiver, per_num, pvalue)
                    ### [2] Downstream targets and TFs
                    max_hop <- object@para$max_hop
                    receptor_tf <- NULL
                    if (nrow(lrdb) > 0) {
                        receptor_tf <- .get_tf_res(celltype_sender, celltype_receiver, lrdb, ggi_tf, cell_pair, st_data, max_hop, co_exp_ratio)
                        if (!is.null(receptor_tf)) {
                            # calculate score
                            lrdb <- .get_score(lrdb, receptor_tf)
                        } else {
                            lrdb <- "NA"
                        }
                    }
                    lrpair <- object@lrpair
                    if (nrow(lrpair) == 0) {
                        if (is.data.frame(lrdb)) {
                            object@lrpair <- lrdb
                        }
                    } else {
                        if (is.data.frame(lrdb)) {
                            lrpair <- rbind(lrpair, lrdb)
                            object@lrpair <- unique(lrpair)
                        }
                    }
                    tf <- object@tf
                    if (nrow(tf) == 0) {
                        if (is.data.frame(receptor_tf)) {
                            object@tf <- receptor_tf
                        }
                    } else {
                        if (is.data.frame(receptor_tf)) {
                            tf <- rbind(tf, receptor_tf)
                            object@tf <- unique(tf)
                        }
                    }
                    object@cellpair[[paste0(celltype_sender, " -- ", celltype_receiver)]] <- cell_pair
                }
            }
            sym <- crayon::combine_styles("bold", "green")
            cat(sym("***Done***"), paste0(celltype_sender, " -- ", celltype_receiver),"\n")
        }
    }
    object@para$min_pairs <- min_pairs
    object@para$min_pairs_ratio <- min_pairs_ratio
    object@para$per_num <- per_num
    object@para$pvalue <- pvalue
    object@para$co_exp_ratio <- co_exp_ratio
    return(object)
}

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
#' @importFrom stats fisher.test
#' @export

get_lr_path <- function(object, celltype_sender, celltype_receiver, ligand, receptor, min_gene_num = 5) {
    # check
    if (!is(object, "SpaTalk")) {
        stop("Invalid class for object: must be 'SpaTalk'!")
    }
    pathways <- object@lr_path$pathways
    ggi_tf <- unique(pathways[, c("src", "dest", "src_tf", "dest_tf")])
    st_data <- .get_st_data(object)
    st_meta <- .get_st_meta(object)
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
    # pathway
    ggi_pathway <- object@lr_path$pathways
    rec_pathway_all <- ggi_pathway[ggi_pathway$src == receptor | ggi_pathway$dest == receptor, ]
    rec_pathway_all <- unique(rec_pathway_all$pathway)
    rec_pathway_yes <- rep("NO", length(rec_pathway_all))
    rec_pathway_gene <- list()
    rec_pathway_pvalue <- list()
    rec_pathway_pvalue <- as.double(rec_pathway_pvalue)
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
            a <- matrix(c(gene_pathway_yes_num, gene_rec_num - gene_pathway_yes_num,
                gene_pathway_num - gene_pathway_yes_num, gene_all_num - gene_rec_num +
                  gene_pathway_yes_num - gene_pathway_num), nrow = 2)
            pvalue <- stats::fisher.test(a)
            pvalue <- as.double(pvalue$p.value)
            rec_pathway_gene[[j]] <- gene_pathway_yes
            rec_pathway_pvalue[j] <- pvalue
            rec_pathway_yes[j] <- "YES"
        }
    }
    rec_pathway_all <- rec_pathway_all[which(rec_pathway_yes == "YES")]
    rec_pathway_pvalue <- rec_pathway_pvalue[which(rec_pathway_yes == "YES")]
    rec_pathway_res <- data.frame(celltype_sender = celltype_sender, celltype_receiver = celltype_receiver,
        ligand = ligand, receptor = receptor, receptor_pathways = rec_pathway_all,
        pvalue = rec_pathway_pvalue, gene_count = 0, stringsAsFactors = FALSE)
    rec_pathway_res$gene <- "NA"
    rec_pathway_gene_yes <- which(rec_pathway_yes == "YES")
    for (i in 1:length(rec_pathway_gene_yes)) {
        genename <- rec_pathway_gene[[rec_pathway_gene_yes[i]]]
        genename1 <- genename[1]
        if (length(genename) > 1) {
            for (j in 2:length(genename)) {
                genename2 <- genename[j]
                genename1 <- paste(genename1, genename2, sep = ",")
            }
        }
        rec_pathway_res$gene_count[i] <- length(genename)
        rec_pathway_res$gene[i] <- genename1
    }
    return(list(tf_path = tf_path_all, path_pvalue = rec_pathway_res))
}
