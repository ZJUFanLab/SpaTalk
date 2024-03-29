% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{dec_cci}
\alias{dec_cci}
\title{Decomposing cell-cell communications for spatial transciptomics data}
\usage{
dec_cci(
  object,
  celltype_sender,
  celltype_receiver,
  n_neighbor = 10,
  min_pairs = 5,
  min_pairs_ratio = 0,
  per_num = 1000,
  pvalue = 0.05,
  co_exp_ratio = 0.1,
  if_doParallel = T,
  use_n_cores = NULL
)
}
\arguments{
\item{object}{SpaTalk object after \code{\link{find_lr_path}}.}

\item{celltype_sender}{Name of celltype_sender.}

\item{celltype_receiver}{Name of celltype_receiver.}

\item{n_neighbor}{Number of neighbor cells to select as the proximal cell-cell pair. Default is \code{10}.}

\item{min_pairs}{Min proximal cell-cell pairs between for sending and receiving cell types. Default is \code{5}.}

\item{min_pairs_ratio}{Min proximal cell-cell pairs ratio between for sending and receiving cell types. Default is \code{0}.}

\item{per_num}{Number of repeat times for permutation test. Default is \code{1000}.}

\item{pvalue}{Include the significantly proximal LR pairs with this cutoff of p value from permutation test. Default is \code{0.05}.}

\item{co_exp_ratio}{Min cell ratio in receiving cells with co-expressed source and target genes for predicting the downstream pathway activity.}

\item{if_doParallel}{Use doParallel. Default is TRUE.}

\item{use_n_cores}{Number of CPU cores to use. Default is all cores - 2.}
}
\value{
SpaTalk object containing the inferred LR pairs and pathways.
}
\description{
Identify the cell-cell communications for single-cell or spot-based spatial transciptomics data with proximal ligand-receptor-target interactions.
}
