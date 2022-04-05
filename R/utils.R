.percent_cell <- function(x) {
    return(length(x[x > 0]))
}

.generate_ref <- function(sc_ndata, sc_celltype) {
    ref_celltype <- unique(sc_celltype$celltype)
    sc_ref <- list()
    for (i in 1:length(ref_celltype)) {
        sc_data_celltype <- sc_ndata[, sc_celltype[sc_celltype$celltype == ref_celltype[i],
            ]$cell]
        if (is(sc_data_celltype, "dgCMatrix")) {
            sc_ref[[i]] <- rowMeans(sc_data_celltype)
            names(sc_ref)[i] <- ref_celltype[i]
        } else {
            names(sc_data_celltype) <- rownames(sc_ndata)
            sc_ref[[i]] <- sc_data_celltype
            names(sc_ref)[i] <- ref_celltype[i]
        }
    }
    sc_ref <- as.data.frame(sc_ref, stringsAsFactors = F)
    cellname <- colnames(sc_ref)
    cellname <- cellname[order(cellname)]
    sc_ref <- sc_ref[, cellname]
    return(sc_ref)
}

.normalize_data <- function(rawdata) {
    rawdata <- Seurat::CreateSeuratObject(rawdata)
    rawdata <- Seurat::NormalizeData(rawdata, verbose = F)
    rawdata <- rawdata[["RNA"]]@data
    return(rawdata)
}

.hvg <- function(sc_ndata, sc_celltype) {
    sc_data1 <- Seurat::CreateSeuratObject(sc_ndata)
    Seurat::Idents(sc_data1) <- sc_celltype$celltype
    sc_hvg <- Seurat::FindAllMarkers(sc_data1)
    sc_hvg <- unique(sc_hvg$gene)
    return(sc_hvg)
}

.coef_nor <- function(st_coef) {
    for (i in 1:nrow(st_coef)) {
        st_coef1 <- as.numeric(st_coef[i, ])
        st_coef1 <- st_coef1/sum(st_coef1)
        st_coef[i, ] <- st_coef1
    }
    return(as.data.frame(st_coef))
}

.st_dist <- function(st_meta) {
    st_dist <- as.matrix(stats::dist(x = cbind(st_meta$x, st_meta$y)))
    rownames(st_dist) <- st_meta[, 1]
    colnames(st_dist) <- st_meta[, 1]
    return(st_dist)
}

.determine_celltype <- function(st_meta, min_percent) {
    st_coef <- st_meta[, -c(1:7)]
    cellname <- colnames(st_coef)
    for (i in 1:nrow(st_coef)) {
        st_coef1 <- as.numeric(st_coef[i, ])
        if (max(st_coef1) > min_percent) {
            st_meta$celltype[i] <- cellname[which(st_coef1 == max(st_coef1))]
        }
    }
    return(st_meta)
}

.generate_newmeta <- function(genename, st_meta, st_dist, min_percent) {
    # generate new data
    newmeta_spot <- NULL
    newmeta_ratio <- NULL
    newmeta_cell <- NULL
    newmeta_x <- NULL
    newmeta_y <- NULL
    st_meta <- st_meta[st_meta$label != "less nFeatures", ]
    cellname <- colnames(st_meta)[-c(1:7)]
    for (i in 1:nrow(st_meta)) {
        spot_name <- st_meta$spot[i]
        spot_x <- st_meta$x[i]
        spot_y <- st_meta$y[i]
        spot_cellnum <- st_meta$cell_num[i]
        spot_percent <- as.numeric(st_meta[i, -c(1:7)])
        spot_percent <- spot_percent * spot_cellnum
        spot_percent_floor <- floor(spot_percent)
        spot_percent_dec <- spot_percent - spot_percent_floor
        # cell_num < 1
        spot_celltype <- which(spot_percent_dec > min_percent)
        k <- 0
        if (length(spot_celltype) > 0) {
            spot_percent <- spot_percent_dec[spot_celltype]
            spot_celltype <- cellname[spot_celltype]
            for (j in 1:length(spot_celltype)) {
                spot_percent1 <- spot_percent[j]
                newmeta_ratio <- c(newmeta_ratio, spot_percent1)
                newmeta_cell <- c(newmeta_cell, spot_celltype[j])
                k <- k + 1
            }
        }
        # cell_num >= 1
        spot_celltype <- which(spot_percent_floor > 0)
        if (length(spot_celltype) > 0) {
            spot_percent <- spot_percent_floor[spot_celltype]
            spot_celltype <- cellname[spot_celltype]
            spot_cell <- NULL
            for (j in 1:length(spot_celltype)) {
                spot_cell <- c(spot_cell, rep(spot_celltype[j], spot_percent[j]))
            }
            for (j in 1:length(spot_cell)) {
                spot_percent1 <- 1
                newmeta_ratio <- c(newmeta_ratio, spot_percent1)
                newmeta_cell <- c(newmeta_cell, spot_cell[j])
                k <- k + 1
            }
        }
        if (k > 0) {
            newmeta_spot <- c(newmeta_spot, rep(spot_name, k))
            if (k == 1) {
                newmeta_x <- c(newmeta_x, spot_x)
                newmeta_y <- c(newmeta_y, spot_y)
            } else {
                st_dist1 <- st_dist[, spot_name]
                st_dist1 <- st_dist1[st_dist1 > 0]
                st_dist1 <- min(st_dist1)
                for (j in 1:k) {
                  st_angle_new <- sample(x = c(0:360), size = 1)
                  st_dist_new <- sample(x = c(0:st_dist1), size = 1)
                  newmeta_x1 <- spot_x + st_dist_new * cos(st_angle_new * pi/180)
                  newmeta_y1 <- spot_y + st_dist_new * sin(st_angle_new * pi/180)
                  newmeta_x <- c(newmeta_x, newmeta_x1)
                  newmeta_y <- c(newmeta_y, newmeta_y1)
                }
            }
        }
    }
    newmeta <- data.frame(spot = newmeta_spot, cell_ratio = newmeta_ratio, celltype = newmeta_cell,
        x = newmeta_x, y = newmeta_y, stringsAsFactors = F)
    return(newmeta)
}

.generate_newmeta_cell <- function(newmeta, st_ndata, sc_ndata, sc_celltype, iter_num) {
    newmeta_spotname <- unique(newmeta$spot)
    newmeta_cell <- NULL
    cat(crayon::cyan("Generating single-cell data for each spot", "\n"))
    pb <- progress::progress_bar$new(format = "[:bar] Finished::percent Time :elapsedfull",
        total = length(newmeta_spotname), clear = FALSE, width = 60, complete = "+",
        incomplete = "-")
    for (i in 1:length(newmeta_spotname)) {
        spot_name <- newmeta_spotname[i]
        newmeta_spot <- .generate_newmeta_spot(spot_name, newmeta, st_ndata, sc_ndata,
            sc_celltype, iter_num)
        newmeta_cell <- rbind(newmeta_cell, newmeta_spot)
        pb$tick()
    }
    newmeta_cell$cell <- paste0("C", 1:nrow(newmeta))
    newmeta_cell <- newmeta_cell[, c(8, 4, 5, 3:1, 7, 6)]
    return(newmeta_cell)
}

.generate_newmeta_spot <- function(spot_name, newmeta, st_ndata, sc_ndata, sc_celltype,
    iter_num) {
    newmeta_spot <- newmeta[newmeta$spot == spot_name, ]
    spot_ndata <- as.numeric(st_ndata[, spot_name])
    # random sampling
    score_cor <- NULL
    spot_cell_id <- list()
    for (k in 1:iter_num) {
        spot_cell_id_k <- NULL
        for (j in 1:nrow(newmeta_spot)) {
            spot_celltype <- newmeta_spot$celltype[j]
            sc_celltype1 <- sc_celltype[sc_celltype$celltype == spot_celltype, "cell"]
            sc_celltype1 <- sc_celltype1[sample(x = 1:length(sc_celltype1), size = 1)]
            spot_cell_id_k <- c(spot_cell_id_k, sc_celltype1)
        }
        if (length(spot_cell_id_k) == 1) {
            spot_ndata_pred <- as.numeric(sc_ndata[, spot_cell_id_k])
        } else {
            spot_ndata_pred <- as.numeric(rowSums(sc_ndata[, spot_cell_id_k]))
        }
        spot_ndata_cor <- cor(spot_ndata, spot_ndata_pred)
        score_cor <- c(score_cor, spot_ndata_cor)
        spot_cell_id[[k]] <- spot_cell_id_k
    }
    spot_cell_id <- spot_cell_id[[which.max(score_cor)]]
    newmeta_spot$cell_id <- spot_cell_id
    newmeta_spot$cor <- max(score_cor)
    return(newmeta_spot)
}

.get_st_meta <- function(object) {
    st_type <- object@para$st_type
    if (st_type == "single-cell") {
        st_meta <- object@meta$rawmeta
        st_meta <- st_meta[st_meta$celltype != "unsure", ]
        st_meta <- st_meta[st_meta$label != "less nFeatures", ]
        if (nrow(st_meta) == 0) {
            stop("No cell types found in rawmeta by excluding the unsure and less nFeatures cells!")
        }
    } else {
        st_meta <- object@meta$newmeta
    }
    return(st_meta)
}

.get_st_data <- function(object) {
    st_type <- object@para$st_type
    if (st_type == "single-cell") {
        st_data <- object@data$rawndata
        st_meta <- object@meta$rawmeta
        st_meta <- st_meta[st_meta$celltype != "unsure", ]
        st_meta <- st_meta[st_meta$label != "less nFeatures", ]
        st_data <- st_data[, st_meta$cell]
        gene_expressed_ratio <- rowSums(as.matrix(st_data))
        st_data <- st_data[which(gene_expressed_ratio > 0), ]
        if (nrow(st_data) == 0) {
            stop("No cell types found in rawmeta by excluding the unsure and less nFeatures cells!")
        }
    } else {
        st_data <- object@data$newdata
        st_meta <- object@meta$newmeta
        gene_expressed_ratio <- rowSums(as.matrix(st_data))
        st_data <- st_data[which(gene_expressed_ratio > 0), ]
        if (nrow(st_data) == 0) {
            stop("No expressed genes in newdata!")
        }
    }
    return(st_data)
}

.get_cellpair <- function(celltype_dist, st_meta, celltype_sender, celltype_receiver,
    n_neighbor) {
    cell_sender <- st_meta[st_meta$celltype == celltype_sender, ]
    cell_receiver <- st_meta[st_meta$celltype == celltype_receiver, ]
    cell_pair <- list()
    for (j in 1:nrow(cell_sender)) {
        celltype_dist1 <- celltype_dist[, cell_sender$cell[j]]
        celltype_dist1 <- celltype_dist1[celltype_dist1 > 0]
        celltype_dist1 <- celltype_dist1[order(celltype_dist1)]
        cell_pair[[j]] <- names(celltype_dist1[1:n_neighbor])
    }
    names(cell_pair) <- cell_sender$cell
    cell_pair <- as.data.frame(cell_pair, stringsAsFactors = F)
    cell_pair <- reshape2::melt(cell_pair, measure.vars = colnames(cell_pair), variable.name = "cell_sender",
        value.name = "cell_receiver")
    cell_pair$cell_sender <- as.character(cell_pair$cell_sender)
    cell_pair <- cell_pair[cell_pair$cell_receiver %in% cell_receiver$cell, ]
    return(cell_pair)
}

.co_exp <- function(x) {
    x_1 <- x[1:(length(x)/2)]
    x_2 <- x[(length(x)/2 + 1):length(x)]
    x_12 <- x_1 * x_2
    x_12_ratio <- length(x_12[x_12 > 0])/length(x_12)
    return(x_12_ratio)
}

.co_exp_p <- function(x) {
    x_real <- x[length(x)]
    x_per <- x[-length(x)]
    x_p <- length(x_per[x_per >= x_real])/length(x_per)
}

.lr_distance <- function(st_data, cell_pair, lrdb, celltype_sender, celltype_receiver,
    per_num, pvalue) {
    pb <- progress::progress_bar$new(format = "[:bar] Finished::percent Time :elapsedfull",
        total = per_num, clear = FALSE, width = 60, complete = "+", incomplete = "-")
    ### [1] LR distance
    lrdb$celltype_sender <- celltype_sender
    lrdb$celltype_receiver <- celltype_receiver
    lrdb$lr_co_exp_num <- 0
    lrdb$lr_co_ratio <- 0
    lrdb$lr_co_ratio_pvalue <- 1
    rownames(lrdb) <- 1:nrow(lrdb)
    # calculate co-expression ratio
    ndata_ligand <- st_data[lrdb$ligand, cell_pair$cell_sender]
    ndata_receptor <- st_data[lrdb$receptor, cell_pair$cell_receiver]
    ndata_lr <- cbind(ndata_ligand, ndata_receptor)
    lrdb$lr_co_ratio <- apply(ndata_lr, 1, .co_exp)
    lrdb$lr_co_exp_num <- apply(ndata_lr, 1, .co_exp) * nrow(cell_pair)
    # permutation test
    res_per <- list()
    for (j in 1:per_num) {
        set.seed(j)
        cell_id <- sample(x = 1:ncol(st_data), size = nrow(cell_pair) * 2, replace = T)
        ndata_ligand <- st_data[lrdb$ligand, cell_id[1:(length(cell_id)/2)]]
        ndata_receptor <- st_data[lrdb$receptor, cell_id[(length(cell_id)/2 + 1):length(cell_id)]]
        ndata_lr <- cbind(ndata_ligand, ndata_receptor)
        res_per[[j]] <- apply(ndata_lr, 1, .co_exp)
        pb$tick()
    }
    names(res_per) <- paste0("P", 1:length(res_per))
    res_per <- as.data.frame(res_per)
    res_per$real <- lrdb$lr_co_ratio
    lrdb$lr_co_ratio_pvalue <- apply(res_per, 1, .co_exp_p)
    lrdb <- lrdb[lrdb$lr_co_ratio_pvalue < pvalue, ]
    return(lrdb)
}

.generate_ggi_res <- function(ggi_tf, cell_pair, receptor_name, st_data, max_hop,
    co_exp_ratio) {
    # generate ggi_res
    ggi_res <- NULL
    ggi_tf1 <- ggi_tf[ggi_tf$src == receptor_name, ]
    ggi_tf1 <- unique(ggi_tf1[ggi_tf1$dest %in% rownames(st_data), ])
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
    # ndata_src and ndata_dest
    ggi_res_temp <- unique(ggi_res[, c("src", "dest")])
    if (nrow(ggi_res_temp) >= 5000) {
        ggi_res <- .co_exp_batch(st_data, ggi_res, cell_pair)
    } else {
        ndata_src <- st_data[ggi_res$src, cell_pair$cell_receiver]
        ndata_dest <- st_data[ggi_res$dest, cell_pair$cell_receiver]
        ndata_gg <- cbind(ndata_src, ndata_dest)
        # calculate co-expression
        ggi_res$co_ratio <- apply(ndata_gg, 1, .co_exp)
    }
    ggi_res <- ggi_res[ggi_res$co_ratio > co_exp_ratio, ]
    return(ggi_res)
}

.co_exp_batch <- function(st_data, ggi_res, cell_pair) {
    ggi_res_temp <- unique(ggi_res[, c("src", "dest")])
    cell_receiver <- unique(cell_pair$cell_receiver)
    m <- floor(nrow(ggi_res_temp)/5000)
    i <- 1
    res <- NULL
    while (i <= (m + 1)) {
        m_int <- 5000 * i
        if (m_int < nrow(ggi_res_temp)) {
            ggi_res_temp1 <- ggi_res_temp[((i - 1) * 5000 + 1):(5000 * i), ]
        } else {
            if (m_int == nrow(ggi_res_temp)) {
                ggi_res_temp1 <- ggi_res_temp[((i - 1) * 5000 + 1):(5000 * i), ]
                i <- i + 1
            } else {
                ggi_res_temp1 <- ggi_res_temp[((i - 1) * 5000 + 1):nrow(ggi_res_temp),
                  ]
            }
        }
        ndata_src <- st_data[ggi_res_temp1$src, cell_receiver]
        ndata_dest <- st_data[ggi_res_temp1$dest, cell_receiver]
        ndata_gg <- cbind(ndata_src, ndata_dest)
        # calculate co-expression
        ggi_res_temp1$co_ratio <- apply(ndata_gg, 1, .co_exp)
        res <- rbind(res, ggi_res_temp1)
        i <- i + 1
    }
    res$merge_key <- paste0(res$src, res$dest)
    ggi_res$merge_key <- paste0(ggi_res$src, ggi_res$dest)
    ggi_res <- merge(ggi_res, res, by = "merge_key", all.x = T, sort = F)
    ggi_res <- ggi_res[, c(2:6, 9)]
    colnames(ggi_res) <- c("src", "dest", "src_tf", "dest_tf", "hop", "co_ratio")
    return(ggi_res)
}

.generate_tf_gene_all <- function(ggi_res, max_hop) {
    tf_gene_all <- NULL
    ggi_hop <- ggi_res[ggi_res$hop == 1, ]
    for (k in 1:max_hop) {
        ggi_hop_yes <- ggi_hop[ggi_hop$dest_tf == "YES", ]
        if (nrow(ggi_hop_yes) > 0) {
            ggi_hop_tf <- ggi_res[ggi_res$hop == k + 1, ]
            if (nrow(ggi_hop_tf) > 0) {
                ggi_hop_yes <- ggi_hop_yes[ggi_hop_yes$dest %in% ggi_hop_tf$src, ]
                if (nrow(ggi_hop_yes) > 0) {
                  tf_gene <- ggi_hop_yes$hop
                  names(tf_gene) <- ggi_hop_yes$dest
                  tf_gene_all <- c(tf_gene_all, tf_gene)
                }
            }
        }
        ggi_hop_no <- ggi_hop[ggi_hop$dest_tf == "NO", ]
        ggi_hop <- ggi_res[ggi_res$hop == k + 1, ]
        ggi_hop <- ggi_hop[ggi_hop$src %in% ggi_hop_no$dest, ]
    }
    return(tf_gene_all)
}

.generate_tf_res <- function(tf_gene_all, celltype_sender, celltype_receiver, receptor_name,
    ggi_res) {
    receptor_tf_temp <- data.frame(celltype_sender = celltype_sender, celltype_receiver = celltype_receiver,
        receptor = receptor_name, tf = names(tf_gene_all), n_hop = as.numeric(tf_gene_all),
        n_target = 0, stringsAsFactors = F)
    tf_names <- names(tf_gene_all)
    tf_n_hop <- as.numeric(tf_gene_all)
    for (i in 1:length(tf_names)) {
        ggi_res_tf <- ggi_res[ggi_res$src == tf_names[i] & ggi_res$hop == tf_n_hop[i] +
            1, ]
        receptor_tf_temp$n_target[i] <- length(unique(ggi_res_tf$dest))
    }
    return(receptor_tf_temp)
}

.get_score <- function(lrdb, receptor_tf) {
    lrdb$score <- 0
    for (j in 1:nrow(lrdb)) {
        receptor_name <- lrdb$receptor[j]
        score_lr <- 1 - lrdb$lr_co_ratio_pvalue[j]
        if (receptor_name %in% receptor_tf$receptor) {
            receptor_tf_temp <- receptor_tf[receptor_tf$receptor == receptor_name,
                ]
            receptor_tf_temp$score_rt <- receptor_tf_temp$n_target * receptor_tf_temp$score/receptor_tf_temp$n_hop
            score_rt <- sum(receptor_tf_temp$score_rt) * (-1)
            score_rt <- 1/(1 + exp(score_rt))
            lrdb$score[j] <- sqrt(score_lr * score_rt)
        }
    }
    lrdb <- lrdb[lrdb$score > 0, ]
    if (nrow(lrdb) == 0) {
        return("NA")
    } else {
        return(lrdb)
    }
}

.get_tf_res <- function(celltype_sender, celltype_receiver, lrdb, ggi_tf, cell_pair,
    st_data, max_hop, co_exp_ratio) {
    receptor_tf <- NULL
    receptor_name <- unique(lrdb$receptor)
    pb <- progress::progress_bar$new(format = "[:bar] Finished::percent Time :elapsedfull",
        total = length(receptor_name), clear = FALSE, width = 60, complete = "+",
        incomplete = "-")
    for (j in 1:length(receptor_name)) {
        # generate ggi_res
        ggi_res <- .generate_ggi_res(ggi_tf, cell_pair, receptor_name[j], st_data,
            max_hop, co_exp_ratio)
        if (nrow(ggi_res) > 0) {
            tf_gene_all <- .generate_tf_gene_all(ggi_res, max_hop)
            tf_gene_all <- data.frame(gene = names(tf_gene_all), hop = tf_gene_all,
                stringsAsFactors = F)
            tf_gene_all_new <- unique(tf_gene_all)
            tf_gene_all <- tf_gene_all_new$hop
            names(tf_gene_all) <- tf_gene_all_new$gene
            ggi_res$dest_tf_enrich <- "NO"
            if (!is.null(tf_gene_all)) {
                ggi_res[ggi_res$dest %in% names(tf_gene_all), ]$dest_tf_enrich <- "YES"
                # generate tf res
                receptor_tf_temp <- .generate_tf_res(tf_gene_all, celltype_sender,
                  celltype_receiver, receptor_name[j], ggi_res)
                # random walk
                receptor_tf_temp$score <- .random_walk(receptor_tf_temp, ggi_res)
                receptor_tf <- rbind(receptor_tf, receptor_tf_temp)
            }
        }
        pb$tick()
    }
    return(receptor_tf)
}

.random_walk <- function(receptor_tf, ggi_res) {
    receptor_name <- unique(receptor_tf$receptor)
    tf_score <- rep(0, nrow(receptor_tf))
    names(tf_score) <- receptor_tf$tf
    max_n_hop <- 10
    for (i in 1:10000) {
        ggi_res1 <- ggi_res[ggi_res$src == receptor_name, ]
        n_hop <- 1
        while (nrow(ggi_res1) > 0 & n_hop <= max_n_hop) {
            target_name <- sample(x = 1:nrow(ggi_res1), size = 1)
            ggi_res1 <- ggi_res1[target_name, ]
            if (ggi_res1$dest %in% names(tf_score)) {
                tf_score[ggi_res1$dest] <- tf_score[ggi_res1$dest] + 1
            }
            ggi_res1 <- ggi_res[ggi_res$src == ggi_res1$dest, ]
            n_hop <- n_hop + 1
        }
    }
    tf_score <- as.numeric(tf_score/10000)
}



.get_tf_path <- function(ggi_res, tf_gene, tf_hop, receptor) {
    tf_path <- NULL
    ggi_res1 <- ggi_res[ggi_res$dest == tf_gene & ggi_res$hop == tf_hop, ]
    if (tf_hop > 1) {
        tf_hop_new <- tf_hop - 1
        for (i in tf_hop_new:1) {
            ggi_res2 <- ggi_res[ggi_res$dest %in% ggi_res1$src & ggi_res$hop == i,
                ]
            ggi_res1 <- ggi_res1[ggi_res1$src %in% ggi_res2$dest, ]
            ggi_res2 <- ggi_res2[ggi_res2$dest %in% ggi_res1$src, ]
            if (i == tf_hop_new) {
                tf_path <- rbind(tf_path, ggi_res1, ggi_res2)
            } else {
                tf_path <- rbind(tf_path, ggi_res2)
            }
            ggi_res1 <- ggi_res2
        }
    } else {
        tf_path <- ggi_res1
    }
    tf_path_new <- NULL
    ggi_res1 <- tf_path[tf_path$src == receptor & tf_path$hop == 1, ]
    if (tf_hop > 1) {
        for (i in 2:tf_hop) {
            ggi_res2 <- tf_path[tf_path$src %in% ggi_res1$dest & tf_path$hop == i,
                ]
            ggi_res2 <- ggi_res2[ggi_res2$src %in% ggi_res1$dest, ]
            if (i == 2) {
                tf_path_new <- rbind(tf_path_new, ggi_res1, ggi_res2)
            } else {
                tf_path_new <- rbind(tf_path_new, ggi_res2)
            }
            ggi_res1 <- ggi_res2
        }
    } else {
        tf_path_new <- ggi_res1
    }
    ggi_res1 <- ggi_res[ggi_res$src == tf_gene & ggi_res$hop == (tf_hop + 1), ]
    tf_path_new <- rbind(tf_path_new, ggi_res1)
    tf_path_new$tf <- tf_gene
    return(tf_path_new)
}
