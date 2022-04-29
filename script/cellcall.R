#==============================================
#================== CellCall ==================
#==============================================
library(tidyverse)
library(cellcall)
options(stringsAsFactors = F)
library(SpaTalk)

load(paste0(system.file(package = 'SpaTalk'),'/extdata/starmap_data.rda'))  # starmap_data
load(paste0(system.file(package = 'SpaTalk'),'/extdata/starmap_meta.rda'))  # starmap_meta

count <- starmap_data
colnames(count) <- str_c(starmap_meta$cell,"|",starmap_meta$celltype)
colnames(count) <- str_replace_all(colnames(count),"_","")
in.content <- count

#1. LOAD DATA
## gene expression stored in the variable in.content
dim(in.content)
in.content[1:4, 1:4]
table(str_split(colnames(in.content), "\\|", simplify = T)[,2])
#2. CREATE OBJECT 
mt <- CreateNichConObject(data=in.content, min.feature = 3,
                          names.field = 2,
                          names.delim = "\\|",
                          source = "CPM",
                          scale.factor = 10^6,
                          Org = species[i],
                          project = "Microenvironment")

mt <- TransCommuProfile(object = mt,
                        pValueCor = 0.05,
                        CorValue = 0.1,
                        topTargetCor=1,
                        p.adjust = 0.05,
                        use.type="median",
                        probs = 0.9,
                        method="weighted",
                        IS_core = TRUE,
                        Org = 'Mus musculus')

res <- mt@data$expr_l_r_log2_scale
res_cellcall <- list()

for (i in colnames(res)) {
  lr <- rownames(res)[res[,i] > 0]
  if(length(lr) != 0){
    ligand <- str_split(lr,"-",simplify = T)[,1]
    receptor <- str_split(lr,"-",simplify = T)[,2]
    res_cellcall[[i]] <- data.frame(ligand = ligand, receptor = receptor, scaled_score = res[res[,i] > 0,i,drop = T])
  }
}

#--- Pathway analysis ---
getActivePATH <- function(mt,receiver,ligand,receptor,Org,IS_core=T){
  if(Org == 'Homo sapiens'){
    f.tmp <- system.file("extdata", "new_ligand_receptor_TFs.txt", package="cellcall")
    triple_df <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
    
    if(IS_core){
    }else{
      f.tmp <- system.file("extdata", "new_ligand_receptor_TFs_extended.txt", package="cellcall")
      triple_relation_extended <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
      triple_df <- rbind(triple_df, triple_relation_extended)
    }
    
    f.tmp <- system.file("extdata", "tf_target.txt", package="cellcall")
    target_df <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
  }else if(Org == 'Mus musculus'){
    f.tmp <- system.file("extdata", "new_ligand_receptor_TFs_homology.txt", package="cellcall")
    triple_df <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
    
    if(IS_core){
    }else{
      f.tmp <- system.file("extdata", "new_ligand_receptor_TFs_homology_extended.txt", package="cellcall")
      triple_relation_extended <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
      triple_df <- rbind(triple_df, triple_relation_extended)
    }
    
    f.tmp <- system.file("extdata", "tf_target_homology.txt", package="cellcall")
    target_df <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
  }
  
  
  
  f.tmp <- system.file("extdata", "KEGG_SYMBOL_ID.txt", package="cellcall")
  pathway.info <- read.table(f.tmp, sep = '\t', quote = "", header = FALSE, stringsAsFactors = F)
  colnames(pathway.info) <- c('id', 'name', "main.object.name", "sub.object.name")
  rownames(pathway.info) <- pathway.info$id
  pathway <- list()
  for(i in 1:length(receptor)){
    active_tf <- names(mt@data$gsea.list[[receiver]]@geneSets)
    triple_tmp <- triple_df %>% filter(Receptor_Symbol == receptor[i] & Ligand_Symbol == ligand[i])
    tf_name <- triple_tmp %>% pull(TF_Symbol)
    if (all(!(tf_name %in% active_tf))) {
      return(0)
    } else {
      tf_name <- tf_name[tf_name %in% active_tf]
    }
    
    target_tmp <- mt@data$gsea.list[[receiver]]@geneSets
    triple_tmp <- triple_tmp %>% filter(TF_Symbol %in% tf_name)
    pathway_id <- triple_tmp %>% pull(pathway_ID)  
    pathway_withoutLR <- data.frame()
    for (k in 1:nrow(triple_tmp)) {
      id_tmp <- unlist(strsplit(pathway_id[k],","))
      if(length(id_tmp) > 1){
        pathway_name <- pathway.info %>% filter(id %in% id_tmp) %>% pull(name)
        pathway_name <- str_c(pathway_name,collapse = ";")
      }else{
        pathway_name <- pathway.info %>% filter(id %in% id_tmp) %>% pull(name)
      }
      r2tf_pair <- data.frame(pathway = pathway_name,source = triple_tmp$Receptor_Symbol[k], target = triple_tmp$TF_Symbol[k])
      tf2t_pair <- data.frame(pathway = rep(pathway_name,length(target_tmp[[triple_tmp$TF_Symbol[k]]])),source = rep(triple_tmp$TF_Symbol[k],length(target_tmp[[triple_tmp$TF_Symbol[k]]])), target = target_tmp[[triple_tmp$TF_Symbol[k]]])
      pathway_withoutLR_tmp <- rbind(r2tf_pair,tf2t_pair)
      pathway_withoutLR <- rbind(pathway_withoutLR,pathway_withoutLR_tmp)
    }
    pathway[[paste(ligand[i],receptor[i],sep = "-")]] <- pathway_withoutLR
  }
  return(pathway)
}

res_cellcall_pathway <- list()
for (i in 1:length(res_cellcall)) {
  
  receiver <- str_split(names(res_cellcall)[i],"-",simplify = T)[1,2]
  ligand <- res_cellcall[[i]]$ligand
  receptor <- res_cellcall[[i]]$receptor
  res_cellcall_pathway[[names(res_cellcall)[i]]] <- getActivePATH(mt,receiver = receiver,ligand = ligand,receptor = receptor,Org = "Mus musculus")
  
}

