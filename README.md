# SpaTalk
### A KG-based cell-cell communication inference tool for spatially resolved transcriptomic data

[![R >4.0](https://img.shields.io/badge/R-%3E%3D4.0-brightgreen)](https://www.r-project.org/) <a href='#devtools'>![installed with devtools](https://img.shields.io/badge/installed%20with-devtools-blue)</a> 

Spatially resolved transcriptomics (ST) provides the informative details of genes and retained the crucial spatial information, which have enabled the uncovering of spatial architecture in intact organs, shedding light on the spatially resolved cell-cell communications mediating tissue homeostasis, development and disease. However, inference of cell-cell communications for ST data remains a great challenge. Here, we present SpaTalk, a spatially resolved cell-cell communication inference method relying on the graph network and knowledge graph to model and score the ligand-receptor-target signaling network between the spatially proximal cells, which were decomposed from the ST data through the non-negative linear model and spatial mapping between single-cell RNA-seq and ST data. SpaTalk is a reliable method that can help scientists uncover the spatially resolved cell-cell communications for either single-cell or spot-based ST data universally, providing new insights into the understanding of spatial cellular dynamics in tissues.

# Workflow
<img src='https://github.com/ZJUFanLab/SpaTalk/blob/main/img/SpaTalk.svg'>

# <a name='devtools'>Install</a>
```
devtools::install_github("ZJUFanLab/SpaTalk")
```
