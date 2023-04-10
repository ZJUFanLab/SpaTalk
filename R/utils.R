.rename_chr <- function(x){
    x <- strsplit(x, split = " ")
    x_new <- NULL
    for (i in 1:length(x)) {
        x1 <- x[[i]]
        if (length(x1) > 1) {
            x2 <- x1[1]
            for (j in 2:length(x1)) {
                x2 <- paste(x2, x1[j], sep = "_")
            }
            x1 <- x2
        }
        x_new <- c(x_new, x1)
    }
    x <- strsplit(x_new, split = "-")
    x_new <- NULL
    for (i in 1:length(x)) {
        x1 <- x[[i]]
        if (length(x1) > 1) {
            x2 <- x1[1]
            for (j in 2:length(x1)) {
                x2 <- paste(x2, x1[j], sep = "_")
            }
            x1 <- x2
        }
        x_new <- c(x_new, x1)
    }
    return(x_new)
}

.show_warning <- function(celltype, celltype_new){
    sc_meta <- data.frame(celltype = celltype, celltype_new = celltype_new, stringsAsFactors = FALSE)
    sc_meta <- unique(sc_meta)
    sc_meta <- sc_meta[sc_meta$celltype != sc_meta$celltype_new, ]
    warning_info <- NULL
    if (nrow(sc_meta) > 0) {
        warning_info <- "celltype of "
        if (nrow(sc_meta) == 1) {
            warning_info <- paste0(warning_info, sc_meta$celltype[1], " has been replaced by ", sc_meta$celltype_new[1])
        } else {
            for (i in 1:nrow(sc_meta)) {
                if (i == nrow(sc_meta)) {
                    warning_info <- paste0(warning_info, "and ",sc_meta$celltype[i])
                } else {
                    warning_info <- paste0(warning_info, sc_meta$celltype[i], ", ")
                }
            }
            warning_info <- paste0(warning_info, " have been replaced by ")
            for (i in 1:nrow(sc_meta)) {
                if (i == nrow(sc_meta)) {
                    warning_info <- paste0(warning_info, "and ",sc_meta$celltype_new[i])
                } else {
                    warning_info <- paste0(warning_info, sc_meta$celltype_new[i], ", ")
                }
            }
        }
    }
    return(warning_info)
}

.percent_cell <- function(x) {
    return(length(x[x > 0]))
}

.run_rctd <- function(st_data, sc_data, st_meta, sc_celltype){
    cellname <- unique(sc_celltype$celltype)
    ref_data <- list()
    ref_celltype <- list()
    for (i in 1:length(cellname)) {
        ref_celltype1 <- sc_celltype[sc_celltype$celltype == cellname[i], ]
        if (nrow(ref_celltype1) < 25) {
            set.seed(i)
            cell_num <- sample(x = 1:nrow(ref_celltype1), size = 25, replace = T)
            cell_new <- ref_celltype1$cell[cell_num]
            ref_data1 <- sc_data[, cell_new]
            ref_celltype1 <- data.frame(cell = colnames(ref_data1), celltype = cellname[i], stringsAsFactors = F)
        } else{
            ref_data1 <- sc_data[, ref_celltype1$cell]
        }
        ref_data[[i]] <- ref_data1
        ref_celltype[[i]] <- ref_celltype1
    }
    ref_data1 <- ref_data[[1]]
    ref_celltype1 <- ref_celltype[[1]]
    for (i in 2:length(cellname)) {
        ref_data1 <- cbind(ref_data1, ref_data[[i]])
        ref_celltype1 <- rbind(ref_celltype1, ref_celltype[[i]])
    }
    ref_celltype1$cell <- paste0("C", 1:nrow(ref_celltype1))
    colnames(ref_data1) <- ref_celltype1$cell
    # reference
    ref_cell_types <- ref_celltype1$celltype
    names(ref_cell_types) <- ref_celltype1$cell
    ref_cell_types <- factor(ref_cell_types)
    ref_nUMI <- colSums(ref_data1)
    ref_refernce <- spacexr::Reference(ref_data1, ref_cell_types, ref_nUMI)
    # test data
    test_spot_coords <- data.frame(xcoord = st_meta$x, ycoord = st_meta$y)
    rownames(test_spot_coords) <- st_meta[,1]
    test_spot_nUMI <- colSums(st_data)
    test_spot_RCTD <- spacexr::SpatialRNA(test_spot_coords, st_data, test_spot_nUMI)
    # run RCTD
    myRCTD <- spacexr::create.RCTD(test_spot_RCTD, ref_refernce)
    myRCTD <- spacexr::run.RCTD(myRCTD, doublet_mode = "full")
    results <- myRCTD@results
    st_coef <- as.matrix(results$weights)
    return(st_coef)
}

.run_seurat <- function(st_data, sc_data, sc_celltype){
    ref_data<- Seurat::CreateSeuratObject(sc_data)
    ref_data<- Seurat::SCTransform(ref_data, verbose = F)
    ref_data<- Seurat::RunPCA(ref_data, verbose = F)
    ref_data$celltype <- sc_celltype$celltype
    # Seurat intergration
    test_spot_data<- Seurat::CreateSeuratObject(st_data)
    test_spot_data<- Seurat::SCTransform(test_spot_data,verbose = F)
    test_spot_data<- Seurat::RunPCA(test_spot_data,verbose = F)
    anchors <- Seurat::FindTransferAnchors(reference = ref_data, query = test_spot_data, normalization.method = "SCT", verbose = F)
    predictions.assay <- Seurat::TransferData(anchorset = anchors, refdata = ref_data$celltype,
        prediction.assay = TRUE, weight.reduction = test_spot_data[["pca"]],dims = 1:30, verbose = F)
    res<- predictions.assay@data
    res<- as.data.frame(t(res))
    res<- res[,-ncol(res)]
    cellname<- colnames(res)
    cellname<- stringr::str_replace_all(cellname,pattern = '-',replacement = '_')
    colnames(res)<- cellname
    st_coef <- as.matrix(res)
    return(st_coef)
}

.run_spotlight <- function(st_data, sc_data, sc_celltype){
    ref_data<- Seurat::CreateSeuratObject(sc_data)
    ref_data<- Seurat::SCTransform(ref_data, verbose = F)
    Seurat::Idents(ref_data)<- sc_celltype$celltype
    cluster_markers_all <- Seurat::FindAllMarkers(object = ref_data, assay = "SCT", slot = "data", verbose = F, only.pos = TRUE)
    ref_data$celltype <- sc_celltype$celltype
    spotlight_ls <- SPOTlight::spotlight_deconvolution(se_sc = ref_data, counts_spatial = as.matrix(st_data), clust_vr = "celltype", cluster_markers = cluster_markers_all)
    spotlight_ls<- as.data.frame(spotlight_ls[[2]])
    spotlight_ls<- spotlight_ls[,-ncol(spotlight_ls)]
    st_coef <- as.matrix(spotlight_ls)
    return(st_coef)
}

.run_deconvSeq <- function(st_data, sc_data, sc_celltype){
    label_sc <- sc_celltype$celltype
    names(label_sc) <- sc_celltype$cell
    label_sc <- as.factor(label_sc)
    design.sc <- model.matrix(~-1+label_sc)
    colnames(design.sc) <- levels(label_sc)
    rownames(design.sc) <- colnames(sc_data)
    dge.sc <- deconvSeq::getdge(as.matrix(sc_data), design.sc, ncpm.min=1, nsamp.min=4, method="bin.loess")
    b0.sc <- deconvSeq::getb0.rnaseq(dge.sc, design.sc, ncpm.min=1, nsamp.min=4)
    dge.st <- deconvSeq::getdge(as.matrix(st_data), NULL, ncpm.min=1, nsamp.min=4, method="bin.loess")
    res <- deconvSeq::getx1.rnaseq(NB0=200, b0.sc, dge.st)
    st_coef <- as.matrix(res$x1)
    return(st_coef)
}

.run_stereoscope <- function(st_data, sc_data, sc_celltype, env, anaconda_path){
    ref_data<- as.data.frame(as.matrix(t(sc_data)))
    ref_data$cell<- rownames(ref_data)
    ref_data<- ref_data[,c(ncol(ref_data),1:(ncol(ref_data)-1))]
    ref_celltype <- sc_celltype
    colnames(ref_celltype)[2] <- 'bio_celltype'
    readr::write_tsv(ref_data,file = 'cnt_data.tsv',quote = "none")
    readr::write_tsv(ref_celltype,file = 'mta_data.tsv',quote = "none")
    test_spot_data<- as.data.frame(as.matrix(t(st_data)))
    test_spot_data$spot<- rownames(test_spot_data)
    test_spot_data<- test_spot_data[,c(ncol(test_spot_data),1:(ncol(test_spot_data)-1))]
    readr::write_tsv(test_spot_data,file = 'cnt_st.tsv',quote = "none")
    if (env == "base") {
      res_sh <- c("#!/bin/sh", paste0(anaconda_path,"/condabin/conda init"), paste0("source ", anaconda_path, "/etc/profile.d/conda.sh"), "conda activate")
        res_sh <- c(res_sh, "stereoscope run --sc_cnt cnt_data.tsv --sc_labels mta_data.tsv -sce 75000  -o res -n 5000 --st_cnt cnt_st.tsv -ste 75000 --gpu -stb 100 -scb 100")
        readr::write_lines(res_sh,file = "res_sh.sh")
    } else {
        res_sh <- c("#!/bin/sh", paste0(anaconda_path,"/condabin/conda init"), paste0("source ", anaconda_path, "/etc/profile.d/conda.sh"), paste0("conda activate ", env))
        res_sh <- c(res_sh, "stereoscope run --sc_cnt cnt_data.tsv --sc_labels mta_data.tsv -sce 75000  -o res -n 5000 --st_cnt cnt_st.tsv -ste 75000 --gpu -stb 100 -scb 100")
        readr::write_lines(x = res_sh,file = "res_sh.sh")
    }
    system("bash res_sh.sh")
    filename <- dir("res/cnt_st/")
    if (length(filename) > 1) {
        for (j in 1:length(filename)) {
            cat(paste0("[", j, "] ", filename[j]),"\n")
        }
        filename_used <- readline(prompt = "Which result is it, select the number: ")
        filename <- filename[as.integer(filename_used)]
    }
    st_coef <- readr::read_tsv(paste0("res/cnt_st/", filename))
    st_coef <- as.data.frame(st_coef)
    rownames(st_coef) <- st_coef[,1]
    st_coef <- st_coef[,-1]
    return(as.matrix(st_coef))
}

.run_cell2location <- function(st_data, st_meta, sc_data, sc_celltype, env){
    sc_meta <- sc_celltype
    rownames(sc_meta) <- sc_meta$cell
    sc_obj <- Seurat::CreateSeuratObject(sc_data,meta.data = sc_meta)
    st_obj <- Seurat::CreateSeuratObject(st_data)
    st_obj[['x']] <- st_meta$x
    st_obj[['y']] <- st_meta$y
    if (!dir.exists("cell2location")) {
        dir.create("cell2location")
    }
    adata_sc <- sceasy::convertFormat(sc_obj, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
    adata_st <- sceasy::convertFormat(st_obj, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
    anndata::write_h5ad(adata_sc,"cell2location/sc_rna.h5ad")
    anndata::write_h5ad(adata_st,"cell2location/st_rna.h5ad")
    res_py <- c("import sys", "import scanpy as sc", "import anndata", "import pandas as pd", "import numpy as np","import matplotlib.pyplot as plt")
    res_py <- c(res_py, "import matplotlib as mpl", "import subprocess", "import time", "import os", "import cell2location", "import scvi")
    res_py <- c(res_py, "from matplotlib import rcParams", "rcParams['pdf.fonttype'] = 42", "import seaborn as sns", "from scipy.sparse import csr_matrix", "from cell2location.utils.filtering import filter_genes")
    res_py <- c(res_py, "sc_file_path = 'cell2location/sc_rna.h5ad'", "spatial_file_path = 'cell2location/st_rna.h5ad'", "celltype_key = 'celltype'", "output_file_path = 'cell2location/'")
    res_py <- c(res_py, "adata_snrna_raw = sc.read_h5ad(sc_file_path)", "adata_vis = sc.read_h5ad(spatial_file_path)", "adata_snrna_raw.X = csr_matrix(adata_snrna_raw.X)", "adata_vis.X = csr_matrix(adata_vis.X)")
    res_py <- c(res_py, "adata_snrna_raw = adata_snrna_raw[~adata_snrna_raw.obs[celltype_key].isin(np.array(adata_snrna_raw.obs[celltype_key].value_counts()[adata_snrna_raw.obs[celltype_key].value_counts() <=1].index))]")
    res_py <- c(res_py, "sc.pp.filter_genes(adata_snrna_raw,min_cells=1)", "sc.pp.filter_cells(adata_snrna_raw,min_genes=1)")
    res_py <- c(res_py, "adata_snrna_raw.obs[celltype_key] = pd.Categorical(adata_snrna_raw.obs[celltype_key])", "adata_snrna_raw = adata_snrna_raw[~adata_snrna_raw.obs[celltype_key].isna(), :]")
    res_py <- c(res_py, "selected = filter_genes(adata_snrna_raw, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)")
    res_py <- c(res_py, "adata_snrna_raw = adata_snrna_raw[:, selected].copy()")
    res_py <- c(res_py, "cell2location.models.RegressionModel.setup_anndata(adata=adata_snrna_raw,labels_key=celltype_key)")
    res_py <- c(res_py, "from cell2location.models import RegressionModel")
    res_py <- c(res_py, "mod = RegressionModel(adata_snrna_raw)")
    res_py <- c(res_py, "mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=True)")
    res_py <- c(res_py, "adata_snrna_raw = mod.export_posterior(adata_snrna_raw, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True})")
    res_py <- c(res_py, "if 'means_per_cluster_mu_fg' in adata_snrna_raw.varm.keys():")
    res_py <- c(res_py, "    inf_aver = adata_snrna_raw.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'", "                                    for i in adata_snrna_raw.uns['mod']['factor_names']]].copy()")
    res_py <- c(res_py, "else:", "    inf_aver = adata_snrna_raw.var[[f'means_per_cluster_mu_fg_{i}'", "                                    for i in adata_snrna_raw.uns['mod']['factor_names']]].copy()")
    res_py <- c(res_py, "inf_aver.columns = adata_snrna_raw.uns['mod']['factor_names']", "inf_aver.iloc[0:5, 0:5]")
    res_py <- c(res_py, "intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)", "adata_vis = adata_vis[:, intersect].copy()", "inf_aver = inf_aver.loc[intersect, :].copy()")
    res_py <- c(res_py, "cell2location.models.Cell2location.setup_anndata(adata = adata_vis)")
    res_py <- c(res_py, "mod = cell2location.models.Cell2location(adata_vis, cell_state_df=inf_aver,N_cells_per_location=30,detection_alpha=200)")
    res_py <- c(res_py, "cell2location.models.Cell2location.view_anndata_setup(mod)", "mod.train(max_epochs=30000,batch_size=None,train_size=1,use_gpu=True)")
    res_py <- c(res_py, "adata_vis = mod.export_posterior(adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True})")
    res_py <- c(res_py, "print(adata_vis)", "adata_vis.obsm['q05_cell_abundance_w_sf'].to_csv(output_file_path + '/cell2location_result.csv')")
    readr::write_lines(x = res_py,file = "res_py.py")
    reticulate::use_condaenv(env)
    reticulate::source_python(file = "res_py.py")
    st_coef <- readr::read_csv("cell2location/cell2location_result.csv")
    st_coef <- as.data.frame(st_coef)
    rownames(st_coef) <- st_coef[,1]
    st_coef <- st_coef[,-1]
    coef_cellname <- colnames(st_coef)
    coef_cellname <- strsplit(coef_cellname,split = "abundance_w_sf_")
    coef_cellname <- as.data.frame(coef_cellname)
    coef_cellname <- as.data.frame(t(coef_cellname))
    coef_cellname <- coef_cellname[,2]
    colnames(st_coef) <- coef_cellname
    return(as.matrix(st_coef))
}

.generate_ref <- function(sc_ndata, sc_celltype, if_doParallel) {
    ref_celltype <- unique(sc_celltype$celltype)
    if (if_doParallel) {
        sc_ref <- foreach::foreach(i = 1:length(ref_celltype), .combine = "cbind", .packages = "Matrix") %dopar% {
            sc_data_celltype <- sc_ndata[, sc_celltype[sc_celltype$celltype == ref_celltype[i], ]$cell]
            if (is(sc_data_celltype, "dgCMatrix")) {
                as.numeric(rowMeans(sc_data_celltype))
            } else {
                as.numeric(sc_data_celltype)
            }
        }
        rownames(sc_ref) <- rownames(sc_ndata)
        colnames(sc_ref) <- ref_celltype
    } else {
        sc_ref <- list()
        for (i in 1:length(ref_celltype)) {
            sc_data_celltype <- sc_ndata[, sc_celltype[sc_celltype$celltype == ref_celltype[i], ]$cell]
            if (is(sc_data_celltype, "dgCMatrix")) {
                sc_ref[[i]] <- rowMeans(sc_data_celltype)
                names(sc_ref)[i] <- ref_celltype[i]
            } else {
                names(sc_data_celltype) <- rownames(sc_ndata)
                sc_ref[[i]] <- sc_data_celltype
                names(sc_ref)[i] <- ref_celltype[i]
            }
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

.generate_newmeta <- function(st_meta, st_dist, min_percent) {
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
                n_neighbor <- 4
                st_dist1 <- st_dist[, spot_name]
                st_dist1 <- st_dist1[st_dist1 > 0]
                st_dist1 <- st_dist1[order(st_dist1)]
                st_dist1 <- st_dist1[1:n_neighbor]
                st_meta_neighbor <- st_meta[st_meta$spot %in% names(st_dist1), ]
                st_meta_neighbor <- .det_neighbor(st_meta_neighbor, spot_x, spot_y, st_dist1)
                if (nrow(st_meta_neighbor) == 0) {
                    for (j in 1:k) {
                        st_angle_new <- sample(x = c(1:360), size = 1)
                        st_dist_new <- sample(x = c(0:st_dist1), size = 1)
                        newmeta_x1 <- spot_x + st_dist_new * cos(st_angle_new * pi/180)
                        newmeta_y1 <- spot_y + st_dist_new * sin(st_angle_new * pi/180)
                        newmeta_x <- c(newmeta_x, newmeta_x1)
                        newmeta_y <- c(newmeta_y, newmeta_y1)
                    }
                } else {
                    for (j in 1:k) {
                        sc_name <- newmeta_cell[j]
                        sc_w1 <- .get_weight1(st_meta_neighbor, sc_name)
                        set.seed(j)
                        st_angle_new <- sample(x = c(1:360), size = 1, prob = as.numeric(sc_w1))
                        spot_ratio <- st_meta[st_meta$spot == spot_name, sc_name]
                        st_dist_new <- .get_weight2(st_meta_neighbor, sc_name, st_dist1, st_angle_new, spot_ratio)
                        newmeta_x1 <- spot_x + st_dist_new * cos(st_angle_new * pi/180)
                        newmeta_y1 <- spot_y + st_dist_new * sin(st_angle_new * pi/180)
                        newmeta_x <- c(newmeta_x, newmeta_x1)
                        newmeta_y <- c(newmeta_y, newmeta_y1)
                    }
                }
            }
        }
    }
    newmeta <- data.frame(spot = newmeta_spot, cell_ratio = newmeta_ratio, celltype = newmeta_cell,
                        x = newmeta_x, y = newmeta_y, stringsAsFactors = F)
    return(newmeta)
}

.get_weight1 <- function(st_meta_neighbor, sc_name){
    sc_name_ratio <- st_meta_neighbor[,sc_name]
    names(sc_name_ratio) <- st_meta_neighbor$w1
    if (!"I" %in% names(sc_name_ratio)) {
        sc_name_ratio <- c(sc_name_ratio, 0)
        names(sc_name_ratio)[length(sc_name_ratio)] <- "I"
    }
    if (!"II" %in% names(sc_name_ratio)) {
        sc_name_ratio <- c(sc_name_ratio, 0)
        names(sc_name_ratio)[length(sc_name_ratio)] <- "II"
    }
    if (!"III" %in% names(sc_name_ratio)) {
        sc_name_ratio <- c(sc_name_ratio, 0)
        names(sc_name_ratio)[length(sc_name_ratio)] <- "III"
    }
    if (!"IV" %in% names(sc_name_ratio)) {
        sc_name_ratio <- c(sc_name_ratio, 0)
        names(sc_name_ratio)[length(sc_name_ratio)] <- "IV"
    }
    sc_name_ratio <- sc_name_ratio[c("I", "II", "III", "IV")]
    sc_name_ratio <- sc_name_ratio + 1
    sc_name_ratio <- sc_name_ratio/sum(sc_name_ratio)
    sc_w1 <- rep(sc_name_ratio, each = 90)
    return(sc_w1)
}

.get_weight2 <- function(st_meta_neighbor, sc_name, st_dist1, st_angle_new, spot_ratio){
    neighbor_weight <- 0
    if (st_angle_new > 0 & st_angle_new <= 90) {
        if ("I" %in% st_meta_neighbor$w1) {
            st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$w1 == "I",]
            neighbor_weight <- st_meta_neighbor1[ ,sc_name]
        }
    }
    if (st_angle_new > 90 & st_angle_new <= 180) {
        if ("II" %in% st_meta_neighbor$w1) {
            st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$w1 == "II",]
            neighbor_weight <- st_meta_neighbor1[ ,sc_name]
        }
    }
    if (st_angle_new > 180 & st_angle_new <= 270) {
        if ("III" %in% st_meta_neighbor$w1) {
            st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$w1 == "III",]
            neighbor_weight <- st_meta_neighbor1[ ,sc_name]
        }
    }
    if (st_angle_new > 270 & st_angle_new <= 360) {
        if ("IV" %in% st_meta_neighbor$w1) {
            st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$w1 == "IV",]
            neighbor_weight <- st_meta_neighbor1[ ,sc_name]
        }
    }
    spot_ratio <- c(spot_ratio, neighbor_weight)
    spot_ratio <- spot_ratio + 1
    spot_ratio <- spot_ratio/sum(spot_ratio)
    sc_w2 <- rep(spot_ratio, each = 5)
    sc_w2 <- sample(c(1:10), size = 1, prob = sc_w2)/10
    st_dist_new <- st_dist1[1]*sc_w2/2
    return(st_dist_new)
}

.det_neighbor <- function(st_meta_neighbor, spot_x, spot_y, st_dist1){
    st_meta_neighbor$w1 <- "NA"
    # right-down
    st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$x >= spot_x & st_meta_neighbor$y > spot_y,]
    if (nrow(st_meta_neighbor1)> 1) {
        neighbor_name <- st_meta_neighbor1$spot
        st_dist2 <- st_dist1[neighbor_name]
        st_dist3 <- names(which(st_dist2 == min(st_dist2)))
        if (length(st_dist3) > 1) {
            st_dist3 <- st_dist3[1]
        }
        st_dist2 <- names(st_dist2)
        st_dist2 <- st_dist2[!st_dist2 %in% st_dist3]
        st_meta_neighbor <- st_meta_neighbor[!st_meta_neighbor$spot %in% st_dist2,]
    }
    st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$x >= spot_x & st_meta_neighbor$y > spot_y,]
    if (nrow(st_meta_neighbor1) > 0) {
        st_meta_neighbor[st_meta_neighbor$spot == st_meta_neighbor1$spot, ]$w1 <- "I"
    }
    # right-down
    st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$x > spot_x & st_meta_neighbor$y <= spot_y,]
    if (nrow(st_meta_neighbor1)> 1) {
        neighbor_name <- st_meta_neighbor1$spot
        st_dist2 <- st_dist1[neighbor_name]
        st_dist3 <- names(which(st_dist2 == min(st_dist2)))
        if (length(st_dist3) > 1) {
            st_dist3 <- st_dist3[1]
        }
        st_dist2 <- names(st_dist2)
        st_dist2 <- st_dist2[!st_dist2 %in% st_dist3]
        st_meta_neighbor <- st_meta_neighbor[!st_meta_neighbor$spot %in% st_dist2,]
    }
    st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$x > spot_x & st_meta_neighbor$y <= spot_y,]
    if (nrow(st_meta_neighbor1) > 0) {
        st_meta_neighbor[st_meta_neighbor$spot == st_meta_neighbor1$spot, ]$w1 <- "IV"
    }
    # left-up
    st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$x < spot_x & st_meta_neighbor$y >= spot_y,]
    if (nrow(st_meta_neighbor1)> 1) {
        neighbor_name <- st_meta_neighbor1$spot
        st_dist2 <- st_dist1[neighbor_name]
        st_dist3 <- names(which(st_dist2 == min(st_dist2)))
        if (length(st_dist3) > 1) {
            st_dist3 <- st_dist3[1]
        }
        st_dist2 <- names(st_dist2)
        st_dist2 <- st_dist2[!st_dist2 %in% st_dist3]
        st_meta_neighbor <- st_meta_neighbor[!st_meta_neighbor$spot %in% st_dist2,]
    }
    st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$x < spot_x & st_meta_neighbor$y >= spot_y,]
    if (nrow(st_meta_neighbor1) > 0) {
        st_meta_neighbor[st_meta_neighbor$spot == st_meta_neighbor1$spot, ]$w1 <- "II"
    }
    # left-down
    st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$x <= spot_x & st_meta_neighbor$y < spot_y,]
    if (nrow(st_meta_neighbor1)> 1) {
        neighbor_name <- st_meta_neighbor1$spot
        st_dist2 <- st_dist1[neighbor_name]
        st_dist3 <- names(which(st_dist2 == min(st_dist2)))
        if (length(st_dist3) > 1) {
            st_dist3 <- st_dist3[1]
        }
        st_dist2 <- names(st_dist2)
        st_dist2 <- st_dist2[!st_dist2 %in% st_dist3]
        st_meta_neighbor <- st_meta_neighbor[!st_meta_neighbor$spot %in% st_dist2,]
    }
    st_meta_neighbor1 <- st_meta_neighbor[st_meta_neighbor$x <= spot_x & st_meta_neighbor$y < spot_y,]
    if (nrow(st_meta_neighbor1) > 0) {
        st_meta_neighbor[st_meta_neighbor$spot == st_meta_neighbor1$spot, ]$w1 <- "III"
    }
    return(st_meta_neighbor)
}

.generate_newmeta_doParallel <- function(st_meta, st_dist, min_percent) {
    # generate new data
    st_meta <- st_meta[st_meta$label != "less nFeatures", ]
    cellname <- colnames(st_meta)[-c(1:7)]
    newmeta <- foreach::foreach (i = 1:nrow(st_meta), .combine = rbind, .packages = "Matrix", .export = c(".det_neighbor", ".get_weight1", ".get_weight2")) %dopar% {
        newmeta_spot <- NULL
        newmeta_ratio <- NULL
        newmeta_cell <- NULL
        newmeta_x <- NULL
        newmeta_y <- NULL
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
                n_neighbor <- 4
                st_dist1 <- st_dist[, spot_name]
                st_dist1 <- st_dist1[st_dist1 > 0]
                st_dist1 <- st_dist1[order(st_dist1)]
                st_dist1 <- st_dist1[1:n_neighbor]
                st_meta_neighbor <- st_meta[st_meta$spot %in% names(st_dist1), ]
                st_meta_neighbor <- .det_neighbor(st_meta_neighbor, spot_x, spot_y, st_dist1)
                if (nrow(st_meta_neighbor) == 0) {
                    for (j in 1:k) {
                        st_angle_new <- sample(x = c(1:360), size = 1)
                        st_dist_new <- sample(x = c(0:st_dist1), size = 1)
                        newmeta_x1 <- spot_x + st_dist_new * cos(st_angle_new * pi/180)
                        newmeta_y1 <- spot_y + st_dist_new * sin(st_angle_new * pi/180)
                        newmeta_x <- c(newmeta_x, newmeta_x1)
                        newmeta_y <- c(newmeta_y, newmeta_y1)
                    }
                } else {
                    for (j in 1:k) {
                        sc_name <- newmeta_cell[j]
                        sc_w1 <- .get_weight1(st_meta_neighbor, sc_name)
                        set.seed(j)
                        st_angle_new <- sample(x = c(1:360), size = 1, prob = as.numeric(sc_w1))
                        spot_ratio <- st_meta[st_meta$spot == spot_name, sc_name]
                        st_dist_new <- .get_weight2(st_meta_neighbor, sc_name, st_dist1, st_angle_new, spot_ratio)
                        newmeta_x1 <- spot_x + st_dist_new * cos(st_angle_new * pi/180)
                        newmeta_y1 <- spot_y + st_dist_new * sin(st_angle_new * pi/180)
                        newmeta_x <- c(newmeta_x, newmeta_x1)
                        newmeta_y <- c(newmeta_y, newmeta_y1)
                    }
                }
            }
            data.frame(spot = newmeta_spot, cell_ratio = newmeta_ratio, celltype = newmeta_cell, x = as.numeric(newmeta_x), y = as.numeric(newmeta_y), stringsAsFactors = F)
        } else {
            data.frame(spot = "NA", cell_ratio = "NA", celltype = "NA", x = "NA", y = "NA", stringsAsFactors = F)
        }
    }
    newmeta <- newmeta[newmeta$spot != "NA",]
    return(newmeta)
}

.generate_newmeta_cell <- function(newmeta, st_ndata, sc_ndata, sc_celltype, iter_num, n_cores, if_doParallel) {
    newmeta_spotname <- unique(newmeta$spot)
    newmeta_cell <- NULL
    cat(crayon::cyan("Generating single-cell data for each spot", "\n"))
    if (if_doParallel) {
        cl <- parallel::makeCluster(n_cores)
        doParallel::registerDoParallel(cl)
        newmeta_cell <- foreach::foreach (i = 1:length(newmeta_spotname), .combine = "rbind", .packages = "Matrix", .export = ".generate_newmeta_spot") %dopar% {
            spot_name <- newmeta_spotname[i]
            .generate_newmeta_spot(spot_name, newmeta, st_ndata, sc_ndata, sc_celltype, iter_num)
        }
        doParallel::stopImplicitCluster()
        parallel::stopCluster(cl)
    } else {
        for (i in 1:length(newmeta_spotname)) {
            spot_name <- newmeta_spotname[i]
            newmeta_spot <- .generate_newmeta_spot(spot_name, newmeta, st_ndata, sc_ndata, sc_celltype, iter_num)
            newmeta_cell <- rbind(newmeta_cell, newmeta_spot)
        }
    }
    newmeta_cell$cell <- paste0("C", 1:nrow(newmeta))
    newmeta_cell <- newmeta_cell[, c(8, 4, 5, 3:1, 7, 6)]
    return(newmeta_cell)
}

.generate_newmeta_spot <- function(spot_name, newmeta, st_ndata, sc_ndata, sc_celltype, iter_num) {
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
    if_skip_dec_celltype <- object@para$if_skip_dec_celltype
    if (st_type == "single-cell") {
        st_meta <- object@meta$rawmeta
        st_meta <- st_meta[st_meta$celltype != "unsure", ]
        st_meta <- st_meta[st_meta$label != "less nFeatures", ]
        if (nrow(st_meta) == 0) {
            stop("No cell types found in rawmeta by excluding the unsure and less nFeatures cells!")
        }
    } else {
        if (if_skip_dec_celltype) {
            st_meta <- object@meta$rawmeta
            colnames(st_meta)[1] <- "cell"
        } else {
            st_meta <- object@meta$newmeta
        }
    }
    return(st_meta)
}

.get_st_data <- function(object) {
    st_type <- object@para$st_type
    if_skip_dec_celltype <- object@para$if_skip_dec_celltype
    if (st_type == "single-cell") {
        st_data <- object@data
        if (if_skip_dec_celltype) {
            st_data <- st_data$rawdata
        } else {
            st_data <- st_data$rawndata
        }
        st_meta <- object@meta$rawmeta
        st_meta <- st_meta[st_meta$celltype != "unsure", ]
        st_meta <- st_meta[st_meta$label != "less nFeatures", ]
        st_data <- st_data[, st_meta$cell]
        gene_expressed_ratio <- rowSums(st_data)
        st_data <- st_data[which(gene_expressed_ratio > 0), ]
        if (nrow(st_data) == 0) {
            stop("No cell types found in rawmeta by excluding the unsure and less nFeatures cells!")
        }
    } else {
        if (if_skip_dec_celltype) {
            st_data <- object@data$rawdata
        } else {
            st_data <- object@data$newdata
        }
        gene_expressed_ratio <- rowSums(st_data)
        st_data <- st_data[which(gene_expressed_ratio > 0), ]
        if (nrow(st_data) == 0) {
            stop("No expressed genes in newdata!")
        }
    }
    return(st_data)
}

.get_cellpair <- function(celltype_dist, st_meta, celltype_sender, celltype_receiver, n_neighbor) {
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
    cell_pair <- reshape2::melt(cell_pair, measure.vars = colnames(cell_pair), variable.name = "cell_sender", value.name = "cell_receiver")
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

.lr_distance_doParallel <- function(st_data, cell_pair, lrdb, celltype_sender, celltype_receiver, per_num, pvalue) {
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
    if (nrow(lrdb) == 1) {
        ndata_lr <- ndata_ligand * ndata_receptor
        lrdb$lr_co_exp_num <- length(ndata_lr[ndata_lr > 0])
        lrdb$lr_co_ratio <- length(ndata_lr[ndata_lr > 0])/length(ndata_lr)
    } else {
        ndata_lr <- cbind(ndata_ligand, ndata_receptor)
        lrdb$lr_co_ratio <- apply(ndata_lr, 1, .co_exp)
        lrdb$lr_co_exp_num <- apply(ndata_lr, 1, .co_exp) * nrow(cell_pair)
    }
    # permutation test
    res_per <- foreach::foreach (j=1:per_num, .packages = "Matrix", .export = ".co_exp") %dopar% {
        set.seed(j)
        cell_id <- sample(x = 1:ncol(st_data), size = nrow(cell_pair) * 2, replace = T)
        ndata_ligand <- st_data[lrdb$ligand, cell_id[1:(length(cell_id)/2)]]
        ndata_receptor <- st_data[lrdb$receptor, cell_id[(length(cell_id)/2 + 1):length(cell_id)]]
        if (nrow(lrdb) == 1) {
            ndata_lr <- ndata_ligand * ndata_receptor
            length(ndata_lr[ndata_lr > 0])/length(ndata_lr)
        } else {
            ndata_lr <- cbind(ndata_ligand, ndata_receptor)
            as.numeric(apply(ndata_lr, 1, .co_exp))
        }
    }
    names(res_per) <- paste0("P", 1:length(res_per))
    res_per <- as.data.frame(res_per)
    res_per$real <- lrdb$lr_co_ratio
    lrdb$lr_co_ratio_pvalue <- apply(res_per, 1, .co_exp_p)
    lrdb <- lrdb[lrdb$lr_co_ratio_pvalue < pvalue, ]
    return(lrdb)
}

.lr_distance <- function(st_data, cell_pair, lrdb, celltype_sender, celltype_receiver, per_num, pvalue) {
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
    ### [1] LR distance
    lrdb$celltype_sender <- celltype_sender
    lrdb$celltype_receiver <- celltype_receiver
    lrdb$lr_co_exp_num <- 0
    lrdb$lr_co_ratio <- 0
    lrdb$lr_co_ratio_pvalue <- 1
    rownames(lrdb) <- 1:nrow(lrdb)
    # calculate co-expression ratio
    if (nrow(lrdb) == 1) {
        ndata_lr <- ndata_ligand * ndata_receptor
        lrdb$lr_co_exp_num <- length(ndata_lr[ndata_lr > 0])
        lrdb$lr_co_ratio <- length(ndata_lr[ndata_lr > 0])/length(ndata_lr)
    } else {
        ndata_ligand <- st_data[lrdb$ligand, cell_pair$cell_sender]
        ndata_receptor <- st_data[lrdb$receptor, cell_pair$cell_receiver]
        ndata_lr <- cbind(ndata_ligand, ndata_receptor)
        lrdb$lr_co_ratio <- apply(ndata_lr, 1, .co_exp)
        lrdb$lr_co_exp_num <- apply(ndata_lr, 1, .co_exp) * nrow(cell_pair)
    }
    # permutation test
    res_per <- list()
    for (j in 1:per_num) {
        set.seed(j)
        cell_id <- sample(x = 1:ncol(st_data), size = nrow(cell_pair) * 2, replace = T)
        ndata_ligand <- st_data[lrdb$ligand, cell_id[1:(length(cell_id)/2)]]
        ndata_receptor <- st_data[lrdb$receptor, cell_id[(length(cell_id)/2 + 1):length(cell_id)]]
        if (nrow(lrdb) == 1) {
            ndata_lr <- ndata_ligand * ndata_receptor
            res_per[[j]] <- length(ndata_lr[ndata_lr > 0])/length(ndata_lr)
        } else {
            ndata_lr <- cbind(ndata_ligand, ndata_receptor)
            res_per[[j]] <- apply(ndata_lr, 1, .co_exp)
        }
    }
    names(res_per) <- paste0("P", 1:length(res_per))
    res_per <- as.data.frame(res_per)
    res_per$real <- lrdb$lr_co_ratio
    lrdb$lr_co_ratio_pvalue <- apply(res_per, 1, .co_exp_p)
    lrdb <- lrdb[lrdb$lr_co_ratio_pvalue < pvalue, ]
    return(lrdb)
}

.generate_ggi_res <- function(ggi_tf, cell_pair, receptor_name, st_data, max_hop, co_exp_ratio) {
    .co_exp <- function(x) {
        x_1 <- x[1:(length(x)/2)]
        x_2 <- x[(length(x)/2 + 1):length(x)]
        x_12 <- x_1 * x_2
        x_12_ratio <- length(x_12[x_12 > 0])/length(x_12)
        return(x_12_ratio)
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
                    ggi_res_temp1 <- ggi_res_temp[((i - 1) * 5000 + 1):nrow(ggi_res_temp), ]
                }
            }
            ndata_src <- st_data[ggi_res_temp1$src, cell_receiver]
            ndata_dest <- st_data[ggi_res_temp1$dest, cell_receiver]
            ndata_gg <- cbind(ndata_src, ndata_dest)
            # calculate co-expression
            ggi_res_temp1$co_ratio <- NA
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
        ggi_res$co_ratio <- NA
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
                ggi_res_temp1 <- ggi_res_temp[((i - 1) * 5000 + 1):nrow(ggi_res_temp), ]
            }
        }
        ndata_src <- st_data[ggi_res_temp1$src, cell_receiver]
        ndata_dest <- st_data[ggi_res_temp1$dest, cell_receiver]
        ndata_gg <- cbind(ndata_src, ndata_dest)
        # calculate co-expression
        ggi_res_temp1$co_ratio <- NA
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

.generate_tf_res <- function(tf_gene_all, celltype_sender, celltype_receiver, receptor_name, ggi_res) {
    receptor_tf_temp <- data.frame(celltype_sender = celltype_sender, celltype_receiver = celltype_receiver,
        receptor = receptor_name, tf = names(tf_gene_all), n_hop = as.numeric(tf_gene_all), n_target = 0, stringsAsFactors = F)
    tf_names <- names(tf_gene_all)
    tf_n_hop <- as.numeric(tf_gene_all)
    for (i in 1:length(tf_names)) {
        ggi_res_tf <- ggi_res[ggi_res$src == tf_names[i] & ggi_res$hop == tf_n_hop[i] + 1, ]
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
            receptor_tf_temp <- receptor_tf[receptor_tf$receptor == receptor_name, ]
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

.get_tf_res_doParallel <- function(celltype_sender, celltype_receiver, lrdb, ggi_tf, cell_pair, st_data, max_hop, co_exp_ratio) {
    receptor_name <- unique(lrdb$receptor)
    receptor_tf <- NULL
    receptor_tf <- foreach::foreach (j=1:length(receptor_name), .packages = "Matrix", .combine = "rbind",
        .export = c(".generate_ggi_res", ".generate_tf_gene_all", ".generate_tf_res", ".random_walk")) %dopar% {
        # generate ggi_res
        ggi_res <- .generate_ggi_res(ggi_tf, cell_pair, receptor_name[j], st_data, max_hop, co_exp_ratio)
        if (nrow(ggi_res) > 0) {
            tf_gene_all <- .generate_tf_gene_all(ggi_res, max_hop)
            tf_gene_all <- data.frame(gene = names(tf_gene_all), hop = tf_gene_all, stringsAsFactors = F)
            tf_gene_all_new <- unique(tf_gene_all)
            tf_gene_all <- tf_gene_all_new$hop
            names(tf_gene_all) <- tf_gene_all_new$gene
            ggi_res$dest_tf_enrich <- "NO"
            if (!is.null(tf_gene_all)) {
                ggi_res[ggi_res$dest %in% names(tf_gene_all), ]$dest_tf_enrich <- "YES"
                # generate tf res
                receptor_tf_temp <- .generate_tf_res(tf_gene_all, celltype_sender, celltype_receiver, receptor_name[j], ggi_res)
                # random walk
                receptor_tf_temp$score <- .random_walk(receptor_tf_temp, ggi_res)
                receptor_tf_temp
            } else {
                data.frame()
            }
        }
        else {
            data.frame()
        }
    }
    return(receptor_tf)
}

.get_tf_res <- function(celltype_sender, celltype_receiver, lrdb, ggi_tf, cell_pair, st_data, max_hop, co_exp_ratio) {
    .co_exp <- function(x) {
        x_1 <- x[1:(length(x)/2)]
        x_2 <- x[(length(x)/2 + 1):length(x)]
        x_12 <- x_1 * x_2
        x_12_ratio <- length(x_12[x_12 > 0])/length(x_12)
        return(x_12_ratio)
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
                    ggi_res_temp1 <- ggi_res_temp[((i - 1) * 5000 + 1):nrow(ggi_res_temp), ]
                }
            }
            ndata_src <- st_data[ggi_res_temp1$src, cell_receiver]
            ndata_dest <- st_data[ggi_res_temp1$dest, cell_receiver]
            ndata_gg <- cbind(ndata_src, ndata_dest)
            # calculate co-expression
            ggi_res_temp1$co_ratio <- NA
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
    .generate_ggi_res <- function(ggi_tf, cell_pair, receptor_name, st_data, max_hop, co_exp_ratio) {
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
            ggi_res$co_ratio <- NA
            ggi_res$co_ratio <- apply(ndata_gg, 1, .co_exp)
        }
        ggi_res <- ggi_res[ggi_res$co_ratio > co_exp_ratio, ]
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
    .generate_tf_res <- function(tf_gene_all, celltype_sender, celltype_receiver, receptor_name, ggi_res) {
        receptor_tf_temp <- data.frame(celltype_sender = celltype_sender, celltype_receiver = celltype_receiver,
                                       receptor = receptor_name, tf = names(tf_gene_all), n_hop = as.numeric(tf_gene_all), n_target = 0, stringsAsFactors = F)
        tf_names <- names(tf_gene_all)
        tf_n_hop <- as.numeric(tf_gene_all)
        for (i in 1:length(tf_names)) {
            ggi_res_tf <- ggi_res[ggi_res$src == tf_names[i] & ggi_res$hop == tf_n_hop[i] + 1, ]
            receptor_tf_temp$n_target[i] <- length(unique(ggi_res_tf$dest))
        }
        return(receptor_tf_temp)
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
    receptor_tf <- NULL
    receptor_name <- unique(lrdb$receptor)
    for (j in 1:length(receptor_name)) {
        # generate ggi_res
        ggi_res <- .generate_ggi_res(ggi_tf, cell_pair, receptor_name[j], st_data, max_hop, co_exp_ratio)
        if (nrow(ggi_res) > 0) {
            tf_gene_all <- .generate_tf_gene_all(ggi_res, max_hop)
            tf_gene_all <- data.frame(gene = names(tf_gene_all), hop = tf_gene_all, stringsAsFactors = F)
            tf_gene_all_new <- unique(tf_gene_all)
            tf_gene_all <- tf_gene_all_new$hop
            names(tf_gene_all) <- tf_gene_all_new$gene
            ggi_res$dest_tf_enrich <- "NO"
            if (!is.null(tf_gene_all)) {
                ggi_res[ggi_res$dest %in% names(tf_gene_all), ]$dest_tf_enrich <- "YES"
                # generate tf res
                receptor_tf_temp <- .generate_tf_res(tf_gene_all, celltype_sender, celltype_receiver, receptor_name[j], ggi_res)
                # random walk
                receptor_tf_temp$score <- .random_walk(receptor_tf_temp, ggi_res)
                receptor_tf <- rbind(receptor_tf, receptor_tf_temp)
            }
        }
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
            ggi_res2 <- ggi_res[ggi_res$dest %in% ggi_res1$src & ggi_res$hop == i, ]
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
            ggi_res2 <- tf_path[tf_path$src %in% ggi_res1$dest & tf_path$hop == i, ]
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

#' @title Show SpaTalk object
#'
#' @return SpaTalk object
#' @import Matrix
#' @importFrom methods show
#'
#' @export

setMethod(
    f = 'show',
    signature = 'SpaTalk',
    definition = function(object) {
        cat("An object of class SpaTalk", "\n")
        st_data <- object@data$rawdata
        st_type <- object@para[["st_type"]]
        lrpair <- object@lrpair
        cat(paste0(nrow(st_data), " genes across ", ncol(st_data), " ", st_type, "s (", nrow(lrpair), " lrpair)"), "\n")
        return(invisible(x = NULL))
    }
)
