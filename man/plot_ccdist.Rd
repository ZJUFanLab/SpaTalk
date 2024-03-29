% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_ccdist}
\alias{plot_ccdist}
\title{Plot cell-cell distribution}
\usage{
plot_ccdist(
  object,
  celltype_sender,
  celltype_receiver,
  color = NULL,
  size = 1,
  if_plot_others = T,
  if_plot_density = T,
  if_plot_edge = T,
  if_show_arrow = T,
  arrow_length = 0.05,
  plot_cells = NULL
)
}
\arguments{
\item{object}{SpaTalk object generated from \code{\link{dec_celltype}}.}

\item{celltype_sender}{Name of celltype_sender.}

\item{celltype_receiver}{Name of celltype_receiver.}

\item{color}{Color for celltype_sender, celltype_receiver, and others. Three values.}

\item{size}{Point size. Default is \code{1}.}

\item{if_plot_others}{Whether to plot others. Default is \code{TRUE}.}

\item{if_plot_density}{Whether to plot marginal density plots. Default is \code{TRUE}.}

\item{if_plot_edge}{Whether to plot edge between neighbors. Default is \code{TRUE}.}

\item{if_show_arrow}{Whether to show the arrow of the plotted edge. Default is \code{TRUE}.}

\item{arrow_length}{Arrow length.}

\item{plot_cells}{Which cells to plot. Default is all cells. Input a character vector of cell names to plot.}
}
\description{
Point plot with spatial distribution of celltype_sender and celltype_receiver
}
