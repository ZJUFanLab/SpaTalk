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

# plot_st_pie()
plot_st_pie(object = obj, pie_scale = 1.3, xy_ratio = 1.8)

plot_st_gene(object = obj,gene = 'Plp1',size = 4, if_use_newmeta = F)
plot_st_gene(object = obj,gene = 'Plp1', if_use_newmeta = T)

plot_st_celltype(object = obj,celltype = 'Oligo')

plot_st_celltype_density(object = obj,celltype = 'Oligo',type = 'raster',color_low = 'purple',color_high = 'yellow')

plot_st_celltype_percent(object = obj,celltype = 'Oligo',size = 4)

plot_st_celltype_all(object = obj)

plot_st_cor_heatmap(object = obj,marker_genes = c("Plp1","Vip","Sst","Lamp5","Pcp4","Mfge8","Pvalb"),
                    celltypes = c("Oligo","VIP","SST","eL2_3","eL6","Astro","PVALB"),scale = "none",
                    if_use_newmeta = F,color_low = 'blue',color_high = 'yellow',color_mid = 'yellow')

# plot_st_cor_heatmap `color_midpoint`
plot_st_cor_heatmap(object = obj,marker_genes = c("Plp1","Vip","Sst"),
                    celltypes = c("eL5","eL6"),scale = "none",if_use_newmeta = T)
