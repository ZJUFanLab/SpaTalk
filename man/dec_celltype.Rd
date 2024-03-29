% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{dec_celltype}
\alias{dec_celltype}
\title{Decomposing cell type for spatial transcriptomics data}
\usage{
dec_celltype(
  object,
  sc_data,
  sc_celltype,
  min_percent = 0.5,
  min_nFeatures = 10,
  if_use_normalize_data = T,
  if_use_hvg = F,
  if_retain_other_genes = F,
  if_doParallel = T,
  use_n_cores = NULL,
  iter_num = 1000,
  method = 1,
  env = "base",
  anaconda_path = "~/anaconda3",
  dec_result = NULL
)
}
\arguments{
\item{object}{SpaTalk object generated from \code{\link{createSpaTalk}}.}

\item{sc_data}{A A data.frame or matrix or dgCMatrix containing counts of single-cell RNA-seq data as the reference, each column representing a cell, each row representing a gene.}

\item{sc_celltype}{A character containing the cell type of the reference single-cell RNA-seq data.}

\item{min_percent}{Min percent to predict new cell type for single-cell st_data or predict new cell for spot-based st_data. Default is \code{0.5}.}

\item{min_nFeatures}{Min number of expressed features/genes for each spot/cell in \code{st_data}. Default is \code{10}.}

\item{if_use_normalize_data}{Whether to use normalized \code{st_data} and \code{sc_data} with Seurat normalization. Default is \code{TRUE}. set it \code{FALSE} when the st_data and sc_data are already normalized matrix with other methods.}

\item{if_use_hvg}{Whether to use highly variable genes for non-negative regression. Default is \code{FALSE}.}

\item{if_retain_other_genes}{Whether to retain other genes which are not overlapped between sc_data and st_data when reconstructing the single-cell ST data. Default is \code{FALSE}. Set it \code{TRUE} to obtain the constructed single-cell ST data with genes consistent with that in sc_data.}

\item{if_doParallel}{Use doParallel. Default is TRUE.}

\item{use_n_cores}{Number of CPU cores to use. Default is all cores - 2.}

\item{iter_num}{Number of iteration to generate the single-cell data for spot-based data. Default is \code{1000}.}

\item{method}{1 means using the SpaTalk deconvolution method, 2 means using RCTD, 3 means using Seurat, 4 means using SPOTlight, 5 means using deconvSeq, 6 means using stereoscope, 7 means using cell2location}

\item{env}{When method set to 6, namely use stereoscope python package to deconvolute, please define the python environment of installed stereoscope. Default is the 'base' environment. Anaconda is recommended. When method set to 7, namely use cell2location python package to deconvolute, please install cell2location to "base" environment.}

\item{anaconda_path}{When using stereoscope, please define the \code{env} parameter as well as the path to anaconda. Default is "~/anaconda3"}

\item{dec_result}{A matrix of deconvolution result from other upcoming methods, row represents spots or cells, column represents cell types of scRNA-seq reference. See \code{\link{demo_dec_result}}}
}
\value{
SpaTalk object containing the decomposing results.
}
\description{
Identify the cellular composition for single-cell or spot-based spatial transcriptomics data with non-negative regression.
}
