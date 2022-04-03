library(SpaTalk)
# As single-cell st_data
load(paste0(system.file(package = 'SpaTalk'), "/extdata/starmap_data.rda"))
load(paste0(system.file(package = 'SpaTalk'), "/extdata/starmap_meta.rda"))

obj <- createSpaTalk(st_data = as.matrix(starmap_data), st_meta = starmap_meta[, -4],
    species = "Mouse", if_st_is_sc = T, spot_max_cell = 1)
obj <- dec_celltype(object = obj, sc_data = as.matrix(starmap_data),sc_celltype = starmap_meta$celltype)
obj@meta$rawmeta$celltype <- starmap_meta$celltype
obj <- find_lr_path(object = obj, lrpairs = lrpairs, pathways = pathways)
obj <- dec_cci_all(object = obj)

# As spot st_data
load(paste0(system.file(package = 'SpaTalk'), "/extdata/starmap_data.rda"))
load(paste0(system.file(package = 'SpaTalk'), "/extdata/starmap_meta.rda"))

st_data <- generate_spot(st_data = as.matrix(starmap_data),st_meta = starmap_meta,
    x_min = 0, x_res = 200, x_max = 3800, y_min = 0, y_res = 500, y_max = 17500)
st_meta <- st_data$st_meta
st_data <- st_data$st_data
obj <- createSpaTalk(st_data = as.matrix(st_data), st_meta = st_meta[, c(1,2,3)],
    species = "Mouse", if_st_is_sc = F, spot_max_cell = 5)
obj <- dec_celltype(object = obj, sc_data = as.matrix(starmap_data),sc_celltype = starmap_meta$celltype)
obj <- find_lr_path(object = obj, lrpairs = lrpairs, pathways = pathways)
obj <- dec_cci_all(object = obj)

# seqFish
devtools::load_all('github_repo/SpaTalk/')

for (i in 1:7) {
    st_meta <- seqfish_ob_meta_raw[seqfish_ob_meta_raw$View == i-1,]
    st_data <- seqfish_ob_data_raw[,st_meta$cell]
}

obj <- createSpaTalk(st_data = as.matrix(st_data), st_meta = st_meta[, c("cell","x","y")],
                     species = "Mouse", if_st_is_sc = T, spot_max_cell = 1)
obj <- dec_celltype(object = obj, sc_data = as.matrix(st_data),sc_celltype = st_meta$celltype)
obj@meta$rawmeta$celltype <- st_meta$celltype
obj <- find_lr_path(object = obj, lrpairs = lrpairs, pathways = pathways)
obj <- dec_cci_all(object = obj)







