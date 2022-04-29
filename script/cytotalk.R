#==============================================
#================== Cytotalk ==================
#==============================================
#cytotalk 4.0.11
library(cytotalk)
library(tidyverse)
options(stringsAsFactors = F)
library(SpaTalk)

# load required dataset
lrp_mouse <-  CytoTalk::lrp_mouse
pcg_mouse <- CytoTalk::pcg_mouse

# reformat & load data

load(paste0(system.file(package = 'SpaTalk'),'/extdata/starmap_data.rda'))  # starmap_data
load(paste0(system.file(package = 'SpaTalk'),'/extdata/starmap_meta.rda'))  # starmap_meta


ndata <- starmap_data
celltype <- starmap_meta
all(colnames(ndata) == celltype$cell)
ndata <- rownames_to_column(ndata,"Gene")
celltype <- celltype %>% select(1,4)
colnames(celltype) <- c("cell","cell_type")
write.table(ndata,file = "cytotalk_analysis/starmap_counts.txt",quote = F,sep = "\t",row.names = F)
write.table(celltype,file = "cytotalk_analysis/starmap_meta.txt",quote = F,sep = "\t",row.names = F)

fpath_mat <- "cytotalk_analysis/starmap_counts.txt"
fpath_meta <- "cytotalk_analysis/starmap_meta.txt"

lst_scrna <- read_matrix_with_meta(fpath_mat, fpath_meta,auto_transform = T)
table(lst_scrna$cell_types)

cellpairs <- combn(unique(lst_scrna$cell_types),2)

# run CytoTalk process
for (i in 1:ncol(cellpairs)) {
  
  type_a <- cellpairs[1,i]
  type_b <- cellpairs[2,i]
  CP_out <- paste(cellpairs[1,i],cellpairs[2,i],sep = "_")
  cat(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>",i,"/",ncol(cellpairs),">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
  cat("run CytoTalk process","\t","This is running for ",CP_out,"\n")
  run_cytotalk(lst_scrna, type_a, type_b,lrp = lrp_mouse, pcg = pcg_mouse,dir_out = paste0("cytotalk_analysis/starmap_result/",CP_out),cores = 10)
  
}

res_cytotalk <- list()
for (i in 1:ncol(cellpairs)){
  cellA <- cellpairs[1,i]
  cellB <- cellpairs[2,i]
  cellAB <- paste(cellA,cellB,sep = "_")
  file.tmp <- "PathwayAnalysis.txt"
  file.abs <- str_c(paste0("cytotalk_analysis/starmap_result/",cellAB,"/"),file.tmp)
  if(file.exists(file.abs)){
    df.tmp <- read.table(file.abs, header = T)
    df.tmp <- df.tmp %>% arrange(pval_potential)
    res_cytotalk[[cellAB]] <- df.tmp
  }
}

all_type <- c()
for (i in 1:ncol(cellpairs)){
cellA <- cellpairs[1,i]
cellB <- cellpairs[2,i]
cellAB <- paste(cellA,cellB,sep = "_")
all_type <- c(all_type,cellAB)
}

res_cytotalk_pathway <- list()
for (i in 1:ncol(cellpairs)){
  cellA <- cellpairs[1,i]
  cellB <- cellpairs[2,i]
  cellAB <- paste(cellA,cellB,sep = "_")
  file.tmp <- dir(paste0("cytotalk_analysis/starmap_result/",cellAB,"/pathways"))
  file.abs <- str_c(paste0("cytotalk_analysis/starmap_result/",cellAB,"/pathways/"),file.tmp)
  if(length(file.tmp) != 0){
    df_lr <- data.frame()
    pathway_list <- list()
    for (k in 1:length(file.tmp)){
      tb_tmp <- read.csv(file.abs[k],sep = "\t")
      tb_tmp <- tb_tmp %>% mutate(node1 = str_extract(node1,"[:alnum:]++"), node2 = str_extract(node2,"[:alnum:]++"))
      lr_tmp <- str_split(file.tmp[k],"\\.")[[1]][1]
      l_tmp <- str_split(lr_tmp,"--",simplify = T)[1]
      r_tmp <- str_split(lr_tmp,"--",simplify = T)[2]
      if(str_detect(l_tmp,"_")){
        sender <- cellB
        receiver <- cellA
      } else if(str_detect(r_tmp,"_")){
        sender <- cellA
        receiver <- cellB
      }
      source <- tb_tmp %>% pull(node1) %>% str_extract("[:alnum:]++")
      target <- tb_tmp %>% pull(node2) %>% str_extract("[:alnum:]++")
      in_cell <- tb_tmp %>% pull(node1_type)
      in_cell <- str_replace(in_cell,sender,paste0("sender_",sender))
      in_cell <- str_replace(in_cell,receiver,paste0("receiver_",receiver))
      index <- which(tb_tmp$is_ct_edge)
      in_cell[index] <- paste0("sender_",sender,"--","receiver_",receiver)
      
      pathway <- file.tmp[k] %>% str_split("\\.",simplify = T) %>% str_replace("_","") %>% .[1]
      df.tmp <- data.frame(pathway = rep(pathway,nrow(tb_tmp)),source = source,target = target,in_cell = in_cell)
      pathway_list[[pathway]] <- df.tmp
    }
    res_cytotalk_pathway[[cellAB]] <- pathway_list
  }
}

# Output res_cytotalk res_cytotalk_pathway