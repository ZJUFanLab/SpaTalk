library(Giotto)
library(tidyverse)
options(stringsAsFactors = F)
library(SpaTalk)


# reformat & load data

load(paste0(system.file(package = 'SpaTalk'),'/extdata/starmap_data.rda'))  # starmap_data
load(paste0(system.file(package = 'SpaTalk'),'/extdata/starmap_meta.rda'))  # starmap_meta

starmap_meta <- starmap_meta %>% mutate(View = 0) %>% rename(c("ID" = "cell","X" = "x", "Y" = "y", "cell_types" = "celltype", "FOV" = "View"))

# 1. set working directory
my_working_dir = 'spatialtalk/starmap'
write.table(starmap_data,'spatialtalk/starmap/starmap_expression.txt')

expr_path <- 'spatialtalk/starmap/starmap_expression.txt'

# 1. (optional) set Giotto instructions
instrs = createGiottoInstructions(save_plot = TRUE, 
                                  show_plot = FALSE,
                                  save_dir = my_working_dir, 
                                  python_path = python_path)

## create file with offset information
my_offset_file = data.table::data.table(field = c(0),
                                        x_offset = c(0),
                                        y_offset = c(0))

## create a stitch file
stitch_file = stitchFieldCoordinates(location_file = celltype,
                                     offset_file = my_offset_file,
                                     cumulate_offset_x = T,
                                     cumulate_offset_y = F,
                                     field_col = 'FOV',
                                     reverse_final_x = F,
                                     reverse_final_y = T)


stitch_file    = stitch_file[,.(ID, X_final, Y_final)]
my_offset_file = my_offset_file[,.(field, x_offset_final, y_offset_final)]

cellmeta <- celltype %>% select(ID,FOV,cell_types)                                     

starmap <- createGiottoObject(raw_exprs = expr_path,
                                 spatial_locs = stitch_file,
                                 offset_file = my_offset_file,
                                 instructions = instrs)

## add additional annotation if wanted
starmap <- addCellMetadata(starmap,
                              new_metadata = cellmeta,
                              by_column = T,
                              column_cell_ID = 'ID')        
                              
## normalize
starmap <- normalizeGiotto(gobject = starmap, scalefactor = 6000, verbose = T)       

## add gene & cell statistics
starmap <- addStatistics(gobject = starmap)

## adjust expression matrix for technical or known variables
starmap <- adjustGiottoMatrix(gobject = starmap, expression_values = c('normalized'),
                                 batch_columns = NULL, covariate_columns = c('nr_genes', 'total_expr'),
                                 return_gobject = TRUE,
                                 update_slot = c('custom'))

starmap <- createSpatialNetwork(gobject = starmap, method = 'kNN', k = 5, name = 'spatial_network')
LR_data = data.table::fread(system.file("extdata", "mouse_ligand_receptors.txt", package = 'Giotto'))

spatial_all_scores <- c()

LR_data[, ligand_det := ifelse(mouseLigand %in% starmap@gene_ID, T, F)]
LR_data[, receptor_det := ifelse(mouseReceptor %in% starmap@gene_ID, T, F)]
LR_data_det = LR_data[ligand_det == T & receptor_det == T]
select_ligands = LR_data_det$mouseLigand
select_receptors = LR_data_det$mouseReceptor

spatial_all_scores = spatCellCellcom(starmap,
                                                spatial_network_name = 'spatial_network',
                                                cluster_column = 'cell_types', 
                                                random_iter = 1000,
                                                gene_set_1 = select_ligands,
                                                gene_set_2 = select_receptors,
                                                adjust_method = 'fdr',
                                                do_parallel = T,
                                                cores = 4,
                                                verbose = 'a little')