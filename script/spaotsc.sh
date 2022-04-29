#-----------------------------------
#--- R script after
#-----------------------------------
# Run before python 
library(SpaTalk)
library(tidyverse)

load(paste0(system.file(package = 'SpaTalk'),'/extdata/starmap_data.rda'))  # starmap_data
load(paste0(system.file(package = 'SpaTalk'),'/extdata/starmap_meta.rda'))  # starmap_meta

ndata <- t(starmap_data)
rownames(ndata) <- NULL

# option for non-view data eg. starmap
starmap_meta <- starmap_meta %>% mutate(View = 0)
cdata.tmp <- starmap_meta %>% dplyr::select(View,cell,x,y)  
cdata.tmp <- cdata.tmp %>% mutate(cell = 1:nrow(cdata.tmp))
cdata.tmp <- cdata.tmp %>% rename(c("cell"='Cell ID',"View"="Field of View","x"="X","y"="Y"))

cluster <- 1:14
names(cluster) <- starmap_meta %>% pull(celltype) %>% unique()
cluster <- cluster[starmap_meta$celltype]
starmap_meta <- starmap_meta %>% mutate(cluster = cluster)

#------ standard operation ------
mdata.tmp <- starmap_meta %>% dplyr::select(cell,cluster)
mdata.tmp <- mdata.tmp %>% mutate(cell = seq_len(nrow(mdata.tmp))-1)
mdata.tmp <- mdata.tmp %>% rename("cell" = "index", "cluster" = "louvain")

write.csv(ndata,"spaotsc/starmap_counts.csv")
write.csv(cdata.tmp,"spaotsc/starmap_cellcentroids.csv",row.names = F)
write.csv(mdata.tmp,"spaotsc/STARMAP_cell_type_annotations.csv",row.names = F)

#===============================================================================================

#-----------------------------------
#--- python script after
#-----------------------------------
from spaotsc import SpaOTsc
import pandas as pd
import numpy as np
from scipy.spatial import distance_matrix
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
import matplotlib as mpl

# please download from https://github.com/zcang/SpaOTsc
ligrec_pairs_lib = np.loadtxt( "/home/lcy/workspace/SpaOTsc/mouse_olfactory_bulb/data/LigandReceptorDatabase/LigRec_secreted_mouse.txt", dtype=str )
def parse_ligrec_pairs(ligrec_pairs):
    ligs = [lig for [lig, rec] in ligrec_pairs]
    recs = [rec for [lig, rec] in ligrec_pairs]

# Read data
datadir_is = "/home/lcy/spatialtalk/spaotsc"
df_count_is = pd.read_csv(datadir_is+"/starmap_counts.csv")
df_centroid_is = pd.read_csv(datadir_is+"/starmap_cellcentroids.csv")

# Get the ligand-receptor pairs measured
genes_measured_is = np.array( df_count_is.columns.values, dtype=str )
ligrec_pairs_is = []
for i in range(len(ligrec_pairs_lib)):
    if ligrec_pairs_lib[i,0] in genes_measured_is and ligrec_pairs_lib[i,1] in genes_measured_is:
        ligrec_pairs_is.append([ligrec_pairs_lib[i,0], ligrec_pairs_lib[i,1]])

# Parse the ligand-receptor information
ligs = []
for [lig, rec] in ligrec_pairs_is:
    if not lig in ligs:
        ligs.append(lig)
recs = {}
for lig in ligs:
    recs_tmp = []
    for i in range(len(ligrec_pairs_is)):
        if ligrec_pairs_is[i][0] == lig:
            recs_tmp.append(ligrec_pairs_is[i][1])
    recs[lig] = recs_tmp
lig_families = ['Adm', 'Agt', 'Angpt', 'Apln', 'Apob', 'Apod', 'Artn', 'Bgn', 'Bmp', 'Hc', 'Ccl', 'Cd', 
                'Cel', 'Chad', 'Clcf', 'Clec', 'Cntn', 'Col', 'Cp', 'Dhh', 'Dkk', 'Ecm', 'Edil', 'Edn',
                'Efn', 'Fbln', 'Fbn', 'Fga', 'Fgf', 'Fn', 'Gas', 'Gdf', 'Gdnf', 'Grem', 'Hhipl', 'Hspg',
                'Ibsp', 'Igf', 'Igfbp', 'Il', 'Inhb', 'Lam', 'Lefty', 'Lep', 'Lgals', 'Lif', 'Lipc', 'Liph',
                'Ltbp', 'Ltf', 'Matn', 'Mmp', 'Mst', 'Mstn', 'Myoc', 'Nampt', 'Ndp', 'Nid', 'Nodal', 'Npnt',
                'Ntn', 'Olfm', 'Osm', 'P4hb', 'Pcsk9', 'Pdgf', 'Pdyn', 'Pla', 'Pltp', 'Pnoc', 'Proc', 'Prss',
                'Qrfp', 'Rbp3', 'Reln', 'Rgma', 'Rspo', 'Sema', 'Serpin', 'Sfrp', 'Slit', 'Sost', 'Tcn', 'Tgf',
                'Thbs', 'Thpo', 'Timp', 'Tnc', 'Tnf', 'Trh', 'Vcan', 'Vegf', 'Vtn', 'Vwf', 'Wnt']
lig_family_members = {}
for lig_family in lig_families:
    lig_list_tmp = []
    for lig in ligs:
        if lig[:len(lig_family)] == lig_family:
            lig_list_tmp.append(lig)
    lig_family_members[lig_family] = lig_list_tmp

# Details taken from the celltype information in metadata.
df_cluster = pd.read_csv(datadir_is+"/STARMAP_cell_type_annotations.csv", sep=',')
louvain_id = np.array( df_cluster["louvain"], int )
#   [1] "eL2_3" "eL6"   "Astro" "PVALB" "Endo"  "VIP"   "SST"   "Smc"   "eL4"   "Micro" "Oligo" "eL5"  
#  [13] "Reln"  "HPC" 
celltype_id = []
for i in louvain_id:
    celltype_id.append(i-1)
celltype_id = np.array(celltype_id)

# Compute the communications
field_index = np.array( df_centroid_is["Field of View"].values, int )
pts_collection = []
S_collection = []
for ifield in range(1):
    df_field = np.log( df_count_is.loc[np.where(field_index==ifield)[0]] + 1 )
    pts = np.array( df_centroid_is[df_centroid_is["Field of View"]==ifield][['X','Y']], float )
    dmat = distance_matrix(pts, pts)
    issc = SpaOTsc.spatial_sc(sc_data=df_field)
    issc.cell_cell_distance(sc_dmat_spatial=dmat)
    #  S = np.zeros_like(dmat)
    S_named = {}
    cellclusterid = celltype_id[np.array(df_centroid_is[df_centroid_is['Field of View'] == ifield].index)]
    n_clusters = len(np.unique(cellclusterid))
    
    for lig_family in lig_families:
        S_cell = np.zeros_like(dmat)
        for lig in lig_family_members[lig_family]:
            for rec in recs[lig]:
                if np.max(issc.sc_data[lig]) < 1.0 or np.sum(issc.sc_data[rec], axis=0) < 1.0: continue  # filter low-expressed LR pair 
                S_cell = issc.spatial_signaling_ot([lig], [rec])
                S_cluster_sample = np.zeros([14,14], float)
                #  calculate score in terms of cluster 
                for clusterA in tuple(np.unique(cellclusterid)):
                    ind_sc_Acl = np.where(cellclusterid==clusterA)[0]
                    for clusterB in tuple(np.unique(cellclusterid)):
                        ind_sc_Bcl = np.where(cellclusterid==clusterB)[0]
                        S_cluster_sample[clusterA,clusterB] = np.mean(S_cell[ind_sc_Acl,:][:,ind_sc_Bcl])
                S_named[lig + '_' + rec] = S_cluster_sample
            #S += S_family
    pts_collection.append(pts)
    S_collection.append(S_named)

# store the flatten S matrix
S_collection_flatten = []

for ifield in range(1):
    S_field = S_collection[ifield]
    S_flatten = {}
    for k,v in S_field.items():
        v = v.flatten()
        S_flatten[k] = v
    S_collection_flatten.append(S_flatten)
for i in range(1):
    csv_to_write = pd.DataFrame.from_dict(S_collection_flatten[i])
    csv_to_write.to_csv("./Field"+str(i)+"_LRscore.csv")

#-----------------------------------
#--- python script before
#-----------------------------------


#=================================================================================


#-----------------------------------
#--- R script after
#-----------------------------------
cell_cluster <- starmap_meta %>% pull(celltype) %>% unique()
cell_cluster_pair <- expand.grid(cell_cluster,cell_cluster)
cell_cluster_pair <- cell_cluster_pair %>% mutate(pairs = paste0(Var2,"|",Var1)) %>% pull(pairs)

LRScore <- list()
for (i in seq_len(1) - 1) {
  LRScore[[i + 1]] <- read.csv(paste0("~/workspace/SpaOTsc/Field", i, "_LRscore.csv"))  # the output file in python
  LRScore[[i + 1]] <- LRScore[[i + 1]][,-1]
  LRScore[[i + 1]] <- t(LRScore[[i + 1]])
  colnames(LRScore[[i + 1]]) <- cell_cluster_pair
  LRScore[[i + 1]] <- as.data.frame(LRScore[[i + 1]])
}

