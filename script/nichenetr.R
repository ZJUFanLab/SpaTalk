#===============================================
#================== NicheNetR ==================
#===============================================
library(SpaTalk)
library(tidyverse)
library(Seurat)
library(nichenetr)

load(paste0(system.file(package = 'SpaTalk'),'/extdata/starmap_data.rda'))  # starmap_data
load(paste0(system.file(package = 'SpaTalk'),'/extdata/starmap_meta.rda'))  # starmap_meta

seuratObj <- CreateSeuratObject(counts = starmap_data)
Idents(seuratObj) <- starmap_meta$celltype

# download NicheNetr database
human_lr_network <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
human_ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
human_weighted_networks <- readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
human_ligand_tf_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_tf_matrix.rds"))
human_sig_network <- readRDS(url("https://zenodo.org/record/3260758/files/signaling_network.rds"))
human_gr_network <- readRDS(url("https://zenodo.org/record/3260758/files/gr_network.rds"))

# We do the Homologous conversion from human gene to mouse gene.
# Then standardize the gene symbol
# Finally we use mouse_xxx_xxx as the variable as replacement of human database

ligand_target_matrix <- mouse_ligand_target_matrix 
weighted_networks <- mouse_weighted_networks 
ligand_tf_matrix <- mouse_ligand_tf_matrix 
sig_network <- mouse_sig_network 
gr_network <- mouse_gr_network
lr_network <- mouse_lr_network

# generate celltype-celltype pair
cellname <- unique(starmap_meta$celltype)

celltype_pair <- NULL
for (i in 1:length(cellname)) {
  d1 <- data.frame(celltype1 = rep(cellname[i], length(cellname)), 
                   celltype2 = cellname,
                   stringsAsFactors = F)
  celltype_pair <- rbind(celltype_pair, d1)
}
celltype_pair <- celltype_pair[celltype_pair$celltype1 != celltype_pair$celltype2, ]


# NicheNet analysis
res <- NULL
res_ggi <- list()
res_ggi_num <- 1
for (i in 1:nrow(celltype_pair)) {
  print(i)
  # sender cell -- mitotic fetal germ cell
  sender_celltypes <- celltype_pair$celltype1[i]
  list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.10)
  expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
  
  # receiver cell -- mural granulosa
  receiver <- celltype_pair$celltype2[i]
  expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.10)
  background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  # DEGs
  DE_table_receiver <- DE_table_receiver_all[DE_table_receiver_all$cluster == receiver,]
  geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
  geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  
  # ligand
  ligands <- lr_network %>% pull(from) %>% unique()
  receptors <- lr_network %>% pull(to) %>% unique()
  
  expressed_ligands <- intersect(ligands, expressed_genes_sender)
  expressed_receptors <- intersect(receptors, expressed_genes_receiver)
  
  potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
  
  # NicheNet ligand activity analysis
  ligand_activities <- predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  ligand_activities <- ligand_activities %>% arrange(-pearson) %>% mutate(rank = 1:nrow(.))
  
  # top 20
  best_upstream_ligands <- ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
  # target
  active_ligand_target_links_df <- best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
  
  active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
  
  order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
  order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
  rownames(active_ligand_target_links) <- rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) <- colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  
  vis_ligand_target <- active_ligand_target_links[order_targets,order_ligands] %>% t()
  
  # Receptors of top-ranked ligands
  lr_network_top <- lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors <- lr_network_top %>% pull(to) %>% unique()
  
  lr_network_top_df_large <- weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  
  lr_network_top_df <- lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix <- lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
  
  dist_receptors <- dist(lr_network_top_matrix, method = "binary")
  hclust_receptors <- hclust(dist_receptors, method = "ward.D2")
  order_receptors <- hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands <- dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands <- hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor <- hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors <- order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor <- order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
  
  vis_ligand_receptor_network <- lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) <- order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) <- order_ligands_receptor %>% make.names()
  
  # strcit
  lr_network_strict <- lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  ligands_bona_fide <- lr_network_strict %>% pull(from) %>% unique()
  receptors_bona_fide <- lr_network_strict %>% pull(to) %>% unique()
  
  lr_network_top_df_large_strict <- lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
  lr_network_top_df_large_strict <- lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))
  
  lr_network_top_df_strict <- lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
  lr_network_top_matrix_strict <- lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)
  
  dist_receptors <- dist(lr_network_top_matrix_strict, method = "binary")
  hclust_receptors <- hclust(dist_receptors, method = "ward.D2")
  order_receptors <- hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands <- dist(lr_network_top_matrix_strict %>% t(), method = "binary")
  hclust_ligands <- hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor <- hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors <- order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
  order_ligands_receptor <- order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))
  
  vis_ligand_receptor_network_strict <- lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network_strict) <- order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network_strict) <- order_ligands_receptor %>% make.names()
  receptors <- rownames(vis_ligand_receptor_network_strict)
  targets <- colnames(vis_ligand_target)
  ligands <- intersect(colnames(vis_ligand_receptor_network_strict),rownames(vis_ligand_target))
  if (length(ligands) > 0) {
    vis_ligand_target<- vis_ligand_target[rownames(vis_ligand_target) %in% ligands,]
    if (length(ligands) == 1) {
      vis_ligand_receptor_network_strict1 <- list()
      vis_ligand_receptor_network_strict1[[1]] <- vis_ligand_receptor_network_strict[,colnames(vis_ligand_receptor_network_strict) == ligands]
      names(vis_ligand_receptor_network_strict1)[1]<- ligands
      vis_ligand_receptor_network_strict <- vis_ligand_receptor_network_strict1
    }
    if (length(ligands) > 1) {
      vis_ligand_receptor_network_strict <- vis_ligand_receptor_network_strict[,colnames(vis_ligand_receptor_network_strict) %in% ligands]
      vis_ligand_receptor_network_strict <- as.list(as.data.frame(vis_ligand_receptor_network_strict))
    }
    lr_pairs<- NULL
    ggi_res<- list()
    for (j in 1:length(vis_ligand_receptor_network_strict)) {
      ligand_name<- names(vis_ligand_receptor_network_strict)[j]
      ligand_data<- vis_ligand_receptor_network_strict[[j]]
      receptor_name<- receptors[which(ligand_data > 0)]
      if (length(vis_ligand_receptor_network_strict) > 1) {
        target_name <- vis_ligand_target[rownames(vis_ligand_target) == ligand_name,]
      } else{
        target_name <- vis_ligand_target
      }
      target_name <- as.numeric(target_name)
      target_name <- targets[which(target_name > 0)]
      if (length(receptor_name) > 0) {
        lr_pair<- data.frame(ligand_gene_symbol = ligand_name,
                             receptor_gene_symbol = receptor_name,stringsAsFactors = F)
        lr_pairs<- rbind(lr_pairs,lr_pair)
      }
      ggi_res[[j]] <- target_name
      names(ggi_res)[j] <- ligand_name
    }
    lr_pairs$celltype_sender <- sender_celltypes
    lr_pairs$celltype_receiver <- receiver
    res <- rbind(res, lr_pairs)
    res_ggi[[res_ggi_num]] <- ggi_res
    names(res_ggi)[res_ggi_num] <- paste(sender_celltypes,'--',receiver)
    res_ggi_num <- res_ggi_num + 1
  }
}

#--- Pathway anaylysis ---
res_pathway <- list()
k <- 0
for(i in names(res_ggi)){
  k <- k + 1 
  print(k)
  for(Ltar in names(res_ggi[[i]])){
    ligands_all <- Ltar
    targets_all <- res_ggi[[i]][[Ltar]]
    targets_all <- str_replace(targets_all,"\\.","-") 
    active_signaling_network <- get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)
    active_signaling_network$sig <- active_signaling_network$sig %>% mutate(type = "sig")
    active_signaling_network$gr <- active_signaling_network$gr %>% mutate(type = "gr")
    active_signaling_network <- rbind(active_signaling_network$sig,active_signaling_network$gr)
    res_pathway[[i]][[ligands_all]] <- active_signaling_network
  }
}

# Output: res res_ggi res_pathway