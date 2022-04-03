# SpaTalk
### A cell-cell communication inference method for spatially resolved transcriptomic data

[![R >4.0](https://img.shields.io/badge/R-%3E%3D4.0-brightgreen)](https://www.r-project.org/) <a href='#devtools'>![installed with devtools](https://img.shields.io/badge/installed%20with-devtools-blue)</a> [![tutorial](https://img.shields.io/badge/tutorial-vignette-yellow)](https://raw.githack.com/ZJUFanLab/SpaTalk/main/vignettes/tutorial.html)

Spatially resolved transcriptomics (ST) provides the informative details of genes and retained the crucial spatial information, which have enabled the uncovering of spatial architecture in intact organs, shedding light on the spatially resolved cell-cell communications mediating tissue homeostasis, development and disease. However, inference of cell-cell communications for ST data remains a great challenge. Here, we present SpaTalk, a spatially resolved cell-cell communication inference method relying on the graph network and knowledge graph to model and score the ligand-receptor-target signaling network between the spatially proximal cells, which were decomposed from the ST data through the non-negative linear model and spatial mapping between single-cell RNA-seq and ST data. SpaTalk is a reliable method that can help scientists uncover the spatially resolved cell-cell communications for either single-cell or spot-based ST data universally, providing new insights into the understanding of spatial cellular dynamics in tissues.

# Workflow
<img src='https://github.com/ZJUFanLab/SpaTalk/blob/main/img/SpaTalk.svg'>

# <a name='devtools'>Install</a>
```
# install devtools and install
install.packages(pkgs = 'devtools')
devtools::install_github('ZJUFanLab/SpaTalk')
```
or
```
# download the source package in Release page and install
# ensure the right directory for SpaTalk-1.0.tar.gz
install.packages(pkgs = 'SpaTalk-1.0.tar.gz',repos = NULL, type = "source")
```
# Usage
SpaTalk method consists of two components, wherein the first is to dissect the cell-type composition of ST data and the second is to infer the spatially resolved cell-cell communications over the decomposed single-cell ST data
- ### Cell-type decomposition to reconstruct single-cell ST atlas with known cell types
```
# object: SpaTalk object containg ST data
# sc_data: A matrix containing counts of scRNA-seq data as the reference
# sc_celltype:  A character containing the cell types for scRNA-seq data

dec_celltype(object, sc_data, sc_celltype)
```

- ### Inference of cell-cell communication and ligand-receptor-target network in space
```
# object: SpaTalk object containg ST and scRNA-seq data
# celltype_sender
# celltype_receiver

dec_cci(object, celltype_sender, celltype_receiver)
```

# Note
[![CellTalkDB v1.0](https://img.shields.io/badge/CellTalkDB-v1.0-blueviolet)](http://tcm.zju.edu.cn/celltalkdb/) [![KEGG pathway](https://img.shields.io/badge/KEGG-pathway-ff69b4)](https://www.kegg.jp/kegg/pathway.html) [![Reactome pathway](https://img.shields.io/badge/Reactome-pathway-brightgreen)](https://reactome.org/) [![AnimalTFDB v3.0](https://img.shields.io/badge/AnimalTFDB-v3.0-yellowgreen)](http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/) 

SpaTalk uses the ligand-receptor interactions (LRIs) from [`CellTalkDB`](http://tcm.zju.edu.cn/celltalkdb/), pathways from [`KEGG`](https://www.kegg.jp/kegg/pathway.html) and [`Reactome`](https://reactome.org/), and transcrptional factors (TFs) from [`AnimalTFDB`](http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/) by default. In the current version:

- __SpaTalk can be applied to either [single-cell (vignette)]() or [spot-based (vignette)]() ST data__
- __SpaTalk allows to use custom [LRIs, pathways, and TFs database (vignette)]()__
- __SpaTalk can visualize [cell-type compositions (vignette)]() and [cell-cell communications (vignette)]()__
- LRIs and pathways can be download at[`data/`](https://github.com/ZJUFanLab/SpaTalk/tree/main/data) 
- Demo data can be download at[`inst/extdata/`](https://github.com/ZJUFanLab/SpaTalk/tree/main/inst/extdata)

__Please refer to the [tutorial vignette](https://raw.githack.com/ZJUFanLab/SpaTalk/main/vignettes/tutorial.html) with demo data processing steps. Detailed functions see the [wiki page](https://github.com/ZJUFanLab/SpaTalk/wiki)__

# About
SpaTalk was developed by Xin Shao. Should you have any questions, please contact Xin Shao at xin_shao@zju.edu.cn
