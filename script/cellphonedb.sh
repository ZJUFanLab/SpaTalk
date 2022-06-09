#-----------------------------------
#--- R script after
#-----------------------------------
# Run before python
library(SpaTalk)
library(tidyverse)

load(paste0(system.file(package = 'SpaTalk'),'/extdata/starmap_data.rda'))  # starmap_data
load(paste0(system.file(package = 'SpaTalk'),'/extdata/starmap_meta.rda'))  # starmap_meta

# orthologs gene
ndata <- starmap_data
orthologs_gene <- mapper[rownames(ndata)]
na_gene <- rownames(ndata)[is.na(orthologs_gene)]
ndata <- ndata[!is.na(orthologs_gene),]
rownames(ndata) <- orthologs_gene[!is.na(orthologs_gene)] %>% make.names()

celltype <- starmap_meta
rownames(celltype) <- celltype$cell
all(colnames(ndata) == celltype$cell)

# create Seurat object 

seuratObj <- CreateSeuratObject(counts = ndata)
Idents(seuratObj) <- celltype$celltype
seuratObj[['celltype']] <- celltype$celltype
seuratObj <- NormalizeData(seuratObj)

writeMM(seuratObj@assays$RNA@counts, file = 'data/cpdb_input/starmap_counts_mtx/matrix.mtx')
# save gene and cell names
write(x = rownames(seuratObj@assays$RNA@counts), file = "data/cpdb_input/starmap_counts_mtx/features.tsv")
write(x = colnames(seuratObj@assays$RNA@counts), file = "data/cpdb_input/starmap_counts_mtx/barcodes.tsv")

seuratObj@meta.data$Cell = rownames(seuratObj@meta.data)
df = seuratObj@meta.data[, c('Cell', 'celltype')]
write.table(df, file ='data/cpdb_input/starmap_meta.tsv', sep = '\t', quote = F, row.names = F)

#===============================================================================================


#-----------------------------------
#--- python script after
#-----------------------------------

conda activate cpdb
cd data/cpdb_input/
cellphonedb method statistical_analysis starmap_meta.tsv starmap_counts_mtx/ --counts-data gene_name

#-----------------------------------
#--- python script before
#-----------------------------------


#=================================================================================

#-----------------------------------
#--- R script after
#-----------------------------------

# standardize cpdb output
res_cpdb_uniondb <- list()
# create mapper to convert gene into mouse gene
mapper <- gene2orthologs$mousegene
names(mapper) <- gene2orthologs$humangene
gene2orthologs <- readRDS("~/待处理/Spatalk_Benchmark/verified_database/cpdb/gene2orthologs_human_mouse.rds")

cpdb2df <- function(data) {
  df_data <- data.frame()
  Ligand = c()
  Receptor = c()
  Ligand.Cluster = c()
  Receptor.Cluster = c()
  isReceptor_fst = c()
  isReceptor_scn = c()
  MeanLR = c()
  
  for (i in 1:nrow(data)) {
    pair <- strsplit(data$`interacting_pair`[i],split = "_")
    for (j in 13:ncol(data)) {
      c_pair = strsplit(colnames(data)[j],"\\|")
      if (data[i,j] != 0.0) {
        Ligand <- c(Ligand,pair[[1]][1])
        Receptor <- c(Receptor,pair[[1]][2])
        Ligand.Cluster <- c(Ligand.Cluster,c_pair[[1]][1])
        Receptor.Cluster <- c(Receptor.Cluster,c_pair[[1]][2])
        isReceptor_fst <- c(isReceptor_fst,data$`receptor_a`[i])
        isReceptor_scn <- c(isReceptor_scn,data$`receptor_b`[i])
        MeanLR <- c(MeanLR,data[i,j])
      }
      
    }
    
  }
  df_data <- data.frame(Ligand = Ligand, Receptor = Receptor,
                        Ligand.Cluster = Ligand.Cluster,Receptor.Cluster = Receptor.Cluster,
                        isReceptor_fst = isReceptor_fst,
                        isReceptor_scn = isReceptor_scn,
                        MeanLR = MeanLR)
  return(df_data)
}

#--- starmap
cpdb_O <- read.table(paste0("cpdb_input/","/out/significant_means.txt"),sep = "\t",header = T,check.names = F,stringsAsFactors = F)
cpdb_O[is.na(cpdb_O)] <- 0
cpdb_O <- cpdb2df(cpdb_O)

cpdb_O <- cpdb_O %>% mutate(Ligand = mapper[Ligand], Receptor = mapper[Receptor])
cpdb_O <- cpdb_O %>% rename(receptor = Ligand,ligand = Receptor, sender = Receptor.Cluster, receiver = Ligand.Cluster)

