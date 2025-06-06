% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/was2code_dist.R
\name{was2code_dist}
\alias{was2code_dist}
\title{Compute Wasserstein-2 Distance Between Individuals}
\usage{
was2code_dist(
  count_input,
  meta_cell,
  meta_ind,
  var_per_cell,
  var2test,
  ncores = 2,
  k = NULL,
  verbose = 0
)
}
\arguments{
\item{count_input}{A count matrix (genes x cells)}

\item{meta_cell}{Cell-level metadata including 'individual' and 'cell_id'}

\item{meta_ind}{Individual-level metadata including 'individual' and test variable}

\item{var_per_cell}{Variable to normalize per cell (e.g., "nCount_RNA")}

\item{var2test}{Binary variable to test at the individual level}

\item{ncores}{Number of cores to use for parallelization}

\item{k}{Number of opposite-group individuals to compare each person with (default = NULL for all)}
}
\value{
A list of 3D distance arrays (one per gene)
}
\description{
Compute Wasserstein-2 Distance Between Individuals
}
