#' @title Demo data of st_data
#'
#' @description Demo data of st_data.
#' @details \code{st_data} used in \code{\link{dec_celltype}} must be a \code{matrix} object, each column representing a spot, each row representing a gene.
#' @return A matrix.
#' @export
#' @examples st_data_demo <- demo_st_data()

demo_st_data <- function() {
    spotname <- paste0("spot", 1:6)
    genename <- c("A1BG", "A2M", "A2MP", "NAT1", "NAT2", "NAT20")
    countdata <- sample(x = c(rep(0, 50), 1:50), size = 36, replace = T)
    st_data <- matrix(countdata, nrow = 6, ncol = 6)
    rownames(st_data) <- genename
    colnames(st_data) <- spotname
    return(st_data)
}

#' @title Demo data of st_meta
#'
#' @description Demo data of st_meta
#' @details \code{st_meta} used in \code{\link{dec_celltype}} must be a \code{data.frame} object with three columns, namely \code{'spot'}, \code{'x'}, \code{'y'} for spot-based spatial transcriptomics data.
#' @return A data.frame.
#' @export
#' @examples st_meta_demo <- demo_st_meta()

demo_st_meta <- function() {
    spotname <- paste0("spot", 1:6)
    x <- sample(1:20, size = 6, replace = T)
    y <- sample(1:20, size = 6, replace = T)
    st_meta <- data.frame(spot = spotname, x = x, y = y, stringsAsFactors = F)
    return(st_meta)
}

#' @title Demo data of single-cell st_data
#'
#' @description Demo data of single-cell st_data.
#' @details \code{st_data} used in \code{\link{dec_celltype}} must be a \code{matrix} object, each column representing a cell, each row representing a gene.
#' @return A matrix.
#' @export
#' @examples st_data_demo <- demo_st_sc_data()

demo_st_sc_data <- function() {
    cellname <- paste0("cell", 1:6)
    genename <- c("A1BG", "A2M", "A2MP", "NAT1", "NAT2", "NAT20")
    countdata <- sample(x = c(rep(0, 50), 1:50), size = 36, replace = T)
    st_data <- matrix(countdata, nrow = 6, ncol = 6)
    rownames(st_data) <- genename
    colnames(st_data) <- cellname
    return(st_data)
}

#' @title Demo data of st_sc_meta
#'
#' @description Demo data of st_sc_meta
#' @details \code{st_sc_meta} used in \code{\link{dec_celltype}} must be a \code{data.frame} object with three columns,  namely \code{'cell'}, \code{'x'}, \code{'y'} for single-cell spatial transcriptomics data.
#' @return A data.frame.
#' @export
#' @examples st_sc_meta_demo <- demo_st_sc_meta()

demo_st_sc_meta <- function() {
    cellname <- paste0("cell", 1:6)
    x <- sample(1:20, size = 6, replace = T)
    y <- sample(1:20, size = 6, replace = T)
    st_sc_meta <- data.frame(cell = cellname, x = x, y = y, stringsAsFactors = F)
    return(st_sc_meta)
}

#' @title Demo data of sc_data
#'
#' @description Demo data of sc_data.
#' @details \code{sc_data} used in \code{\link{dec_celltype}} must be a \code{matrix} object, each column representing a cell, each row representing a gene.
#' @return A matrix.
#' @export
#' @examples sc_data_demo <- demo_sc_data()

demo_sc_data <- function() {
    cellname <- paste0("cell", 1:6)
    genename <- c("A1BG", "A2M", "A2MP", "NAT1", "NAT2", "NATP")
    countdata <- sample(x = c(rep(0, 50), 1:50), size = 36, replace = T)
    sc_data <- matrix(countdata, nrow = 6, ncol = 6)
    rownames(sc_data) <- genename
    colnames(sc_data) <- cellname
    return(sc_data)
}

#' @title Demo data of geneinfo
#'
#' @description Demo data of geneinfo
#' @details \code{geneinfo} used in \code{\link{dec_celltype}} must be a \code{data.frame} object with three columns, namely \code{'symbol'}, \code{'synonyms'}, \code{'species'}.
#' @export
#' @examples geneinfo_demo <- demo_geneinfo()

demo_geneinfo <- function() {
    gene1 <- c("A1BG", "A1BG", "A2MP1", "Aco1")
    gene2 <- c("A1B", "ABG", "A2MP", "Aco")
    species <- c("Human", "Human", "Human", "Mouse")
    geneinfo_demo <- data.frame(symbol = gene1, synonyms = gene2, species = species,
        stringsAsFactors = F)
    return(geneinfo_demo)
}

#' @title Demo data of lrpairs
#'
#' @description Demo data of lrpairs
#' @details \code{lrpairs} used in \code{\link{dec_cci}} must be a \code{data.frame} object with three columns, namely \code{'ligand'}, \code{'receptor'}, \code{'species'}.
#' @return A data.frame.
#' @export
#' @examples lrpairs_demo <- demo_lrpairs()

demo_lrpairs <- function() {
    ligand <- c("CX3CL1", "TGFB1", "CCL2", "Sst")
    receptor <- c("CX3CR1", "TGFBR2", "CCR2", "Sstr1")
    species <- c("Human", "Human", "Human", "Mouse")
    lrpairs_demo <- data.frame(ligand = ligand, receptor = receptor, species = species,
        stringsAsFactors = F)
    return(lrpairs_demo)
}

#' @title Demo data of pathways
#'
#' @description Demo data of pathways
#' @details \code{pathways} used in \code{\link{dec_cci}} must be a \code{data.frame} object with seven columns, namely \code{'src'}, \code{'dest'}, \code{'pathway'}, \code{'source'}, \code{'type'}, \code{'src_tf'}, \code{'dest_tf'}, \code{'species'}.
#' @return A data.frame.
#' @export
#' @examples pathways_demo <- demo_pathways()

demo_pathways <- function() {
    src <- c("CDKN1A", "CDKN1A", "CDK2", "Akt1")
    dest <- c("CDK2", "CDK4", "TP53", "Atf2")
    pathway <- c("p53 signaling pathway", "p53 signaling pathway", "p53 signaling pathway",
        "PI3K-Akt signaling pathway")
    type <- c("Process(activation)", "Process(activation)", "Process(binding)",
        "Process(association)")
    sourcename <- rep("KEGG", 4)
    src_tf <- c("NO", "NO", "NO", "NO")
    dest_tf <- c("NO", "NO", "YES", "YES")
    species <- c("Human", "Human", "Human", "Mouse")
    pathways_demo <- data.frame(src = src, dest = dest, pathway = pathway, source = sourcename,
        type = type, src_tf = src_tf, dest_tf = dest_tf, species = species, stringsAsFactors = F)
    return(pathways_demo)
}

#' @title Demo data of dec_result
#'
#' @description Demo data of dec_result
#' @details \code{dec_result} used in \code{\link{dec_celltype}} must be a \code{matrix} object, each row representing a spot, each column representing a cell type.
#' @return A matrix.
#' @export
#' @examples dec_result_demo <- demo_dec_result()

demo_dec_result <- function(){
    res_weight <- sample(1:100, 18, F)/100
    dec_result_demo <- matrix(res_weight,nrow = 3,ncol = 6)
    colnames(dec_result_demo) <- c("B", "T", "NK", "Monocyte", "Neutrophil", "DC")
    rownames(dec_result_demo) <- paste0("spot", 1:3)
    return(dec_result_demo)
}
