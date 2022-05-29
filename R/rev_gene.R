#' @title Pre-processing step: revising gene symbols
#'
#' @description Revise genes according to NCBI Gene symbols updated in June 30, 2021 for count matrix, user-custom lrpairs data.frame, and user-custom pathways data.frame.
#' @param data A data.frame or matrix or dgCMatrix containing count data each column representing a spot or a cell, each row representing a gene; Or a data.frame containing ligand-receptor pairs; Or a data.frame containing gene-gene interactions and pathways from KEGG and Reactome as well as the information of transcriptional factors.
#' @param data_type A character to define the type of \code{data}, select \code{'count'} for the data matrix, \code{'lrpairs'} for the data.frame containing lrpairs, \code{'pathways'} for the data.frame containing pathways.
#' @param species Species of the data.\code{'Human'} or \code{'Mouse'}.
#' @param geneinfo A data.frame of the system data containing gene symbols of \code{'Human'} and \code{'Mouse'} updated on June 30, 2021 for revising count matrix.
#' @return A new matrix or data.frame.
#' @importFrom crayon red cyan green
#' @import Matrix
#' @export

rev_gene <- function(data = NULL, data_type = NULL, species = NULL, geneinfo = NULL) {
    if (is.null(data)) {
        stop("Please provide the data for revsing gene symbols!")
    }
    if (is.null(data_type) | !is.character(data_type)) {
        stop("Please provide a correct data_type, i.e., 'count', 'lrpairs', or 'pathways'!")
    }
    if (is.null(geneinfo)) {
        stop("Please provide geneinfo for revsing gene symbols, or use the system data like 'geneinfo = geneinfo'")
    }
    if (length(data_type) > 1 | !data_type %in% c("count", "lrpairs", "pathways")) {
        stop("Please provide a correct data_type, i.e., 'count', 'lrpairs', or 'pathways'!")
    }
    if (length(species) > 1 | !species %in% c("Human", "Mouse")) {
        stop("Please provide a correct species, i.e., 'Human', or 'Mouse'!")
    }
    # define the species
    if (species == "Human") {
        geneinfo <- geneinfo[geneinfo$species == "Human", ]
    }
    if (species == "Mouse") {
        geneinfo <- geneinfo[geneinfo$species == "Mouse", ]
    }
    # revise matrix
    if (data_type == "count") {
        if (is(data, "data.frame")) {
            data <- methods::as(as.matrix(data), "dgCMatrix")
        }
        if (is(data, "matrix")) {
            data <- methods::as(data, "dgCMatrix")
        }
        if (!is(data, "dgCMatrix")) {
            stop("st_data must be a data.frame or matrix or dgCMatrix!")
        }
        cat(crayon::cyan("Revising gene symbols for count data", "\n"))
        Sys.sleep(1)
        # revise gene symbols
        genename <- rownames(data)
        genename1 <- genename[genename %in% geneinfo$symbol]
        if (length(genename1) == 0) {
            stop("Please ensure the rownames of data are gene symbols! See demo_st_data()!")
        }
        genename2 <- genename[!genename %in% geneinfo$symbol]
        if (length(genename2) > 0) {
            genename3 <- genename2[genename2 %in% geneinfo$synonyms]
            if (length(genename3) > 0) {
                genename4 <- rep("NA", length(genename3))
                for (i in 1:length(genename3)) {
                  d1 <- geneinfo[geneinfo$synonyms == genename3[i], ]$symbol
                  if (length(d1) == 1) {
                    genename4[i] <- d1
                  }
                }
                genename3 <- c(genename1, genename3)
                genename4 <- c(genename1, genename4)
                genedata <- data.frame(raw_name = genename3, new_name = genename4,
                  stringsAsFactors = F)
                genedata <- genedata[!genedata$new_name == "NA", ]
                genedata1 <- as.data.frame(table(genedata$new_name), stringsAsFactors = F)
                genedata1 <- genedata1[genedata1$Freq == 1, ]
                genedata <- genedata[genedata$new_name %in% genedata1$Var1, ]
                data <- data[genedata$raw_name, ]
                rownames(data) <- genedata$new_name
            }
        } else {
            data <- data[rownames(data) %in% geneinfo$symbol, ]
        }
        cat(crayon::green("***Done***", "\n"))
        Sys.sleep(2)
    }
    # revise lrpairs
    if (data_type == "lrpairs") {
        if (!is.data.frame(data)) {
            stop("data must be a data.frame when data_type is 'lrpairs'!")
        }
        cat(crayon::cyan("Revising gene symbols for lrpairs data.frame", "\n"))
        Sys.sleep(1)
        if (all(c("ligand", "receptor", "species") %in% colnames(data))) {
            # ligand
            genename <- unique(data$ligand)
            genename1 <- genename[genename %in% geneinfo$symbol]
            if (length(genename1) == 0) {
                stop("Please ensure the ligand of data are gene symbols! See demo_lrpairs()!")
            }
            genename2 <- genename[!genename %in% geneinfo$symbol]
            genename2 <- genename[!genename %in% geneinfo$symbol]
            if (length(genename2) > 0) {
                genename3 <- genename2[genename2 %in% geneinfo$synonyms]
                if (length(genename3) > 0) {
                  for (i in 1:length(genename3)) {
                    d1 <- geneinfo[geneinfo$synonyms == genename3[i], ]$symbol
                    if (length(d1) == 1) {
                      data[data$ligand == genename3[i], ]$ligand <- d1
                    }
                  }
                }
            }
            data <- data[data$ligand %in% geneinfo$symbol, ]
            # receptor
            genename <- unique(data$receptor)
            genename1 <- genename[genename %in% geneinfo$symbol]
            if (length(genename1) == 0) {
                stop("Please ensure the receptor of data are gene symbols! See demo_lrpairs()!")
            }
            genename2 <- genename[!genename %in% geneinfo$symbol]
            genename2 <- genename[!genename %in% geneinfo$symbol]
            if (length(genename2) > 0) {
                genename3 <- genename2[genename2 %in% geneinfo$synonyms]
                if (length(genename3) > 0) {
                  for (i in 1:length(genename3)) {
                    d1 <- geneinfo[geneinfo$synonyms == genename3[i], ]$symbol
                    if (length(d1) == 1) {
                      data[data$receptor == genename3[i], ]$receptor <- d1
                    }
                  }
                }
            }
            data <- data[data$receptor %in% geneinfo$symbol, ]
        } else {
            stop("Please provide a correct lrpairs data.frame! See demo_lrpairs()!")
        }
        cat(crayon::green("***Done***", "\n"))
        Sys.sleep(2)
    }
    # revise pathways
    if (data_type == "pathways") {
        if (!is.data.frame(data)) {
            stop("data must be a data.frame when data_type is 'pathways'!")
        }
        cat(crayon::cyan("Revising gene symbols for pathways data.frame", "\n"))
        Sys.sleep(1)
        if (all(c("src", "dest", "pathway", "type", "src_tf", "dest_tf", "species") %in%
            colnames(data))) {
            # src
            genename <- unique(data$src)
            genename1 <- genename[genename %in% geneinfo$symbol]
            if (length(genename1) == 0) {
                stop("Please ensure the src of data are gene symbols! See demo_lrpairs()!")
            }
            genename2 <- genename[!genename %in% geneinfo$symbol]
            genename2 <- genename[!genename %in% geneinfo$symbol]
            if (length(genename2) > 0) {
                genename3 <- genename2[genename2 %in% geneinfo$synonyms]
                if (length(genename3) > 0) {
                  for (i in 1:length(genename3)) {
                    d1 <- geneinfo[geneinfo$synonyms == genename3[i], ]$symbol
                    if (length(d1) == 1) {
                      data[data$src == genename3[i], ]$src <- d1
                    }
                  }
                }
            }
            data <- data[data$src %in% geneinfo$symbol, ]
            # dest
            genename <- unique(data$dest)
            genename1 <- genename[genename %in% geneinfo$symbol]
            if (length(genename1) == 0) {
                stop("Please ensure the dest of data are gene symbols! See demo_pathways()!")
            }
            genename2 <- genename[!genename %in% geneinfo$symbol]
            genename2 <- genename[!genename %in% geneinfo$symbol]
            if (length(genename2) > 0) {
                genename3 <- genename2[genename2 %in% geneinfo$synonyms]
                if (length(genename3) > 0) {
                  for (i in 1:length(genename3)) {
                    d1 <- geneinfo[geneinfo$synonyms == genename3[i], ]$symbol
                    if (length(d1) == 1) {
                      data[data$dest == genename3[i], ]$dest <- d1
                    }
                  }
                }
            }
            data <- data[data$dest %in% geneinfo$symbol, ]
        } else {
            stop("Please provide a correct pathways data.frame! See demo_pathways()!")
        }
        cat(crayon::green("***Done***", "\n"))
        Sys.sleep(2)
    }
    return(data)
}
