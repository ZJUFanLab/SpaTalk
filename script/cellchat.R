#==============================================
#================== CellChat ==================
#==============================================
library(CellChat)
library(SpaTalk)
library(tidyverse)
library(Seurat)
options(stringsAsFactors = FALSE)

load(paste0(system.file(package = 'SpaTalk'),'/extdata/starmap_data.rda'))  # starmap_data
load(paste0(system.file(package = 'SpaTalk'),'/extdata/starmap_meta.rda'))  # starmap_meta
  

ndata <- starmap_data
celltype <- starmap_meta
rownames(celltype) <- celltype$cell
all(colnames(ndata) == celltype$cell)

seuratObj <- CreateSeuratObject(counts = ndata)
seuratObj <- NormalizeData(seuratObj)
Idents(seuratObj) <- celltype$celltype
  
data.input <- GetAssayData(seuratObj, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(seuratObj)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
  
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
CellChatDB <- CellChat::CellChatDB.mouse
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- res_cellchat[[use_db]][['starmap']][['starmap']] <- df.net %>% dplyr::select(1:6)
