% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fd_immigration.R
\name{fd_immigration}
\alias{fd_immigration}
\title{Frequency-dependent immigration}
\usage{
fd_immigration(neighborhood_df, disp_inds, delta, discrete = FALSE, no_spp = 2)
}
\arguments{
\item{neighborhood_df}{neighborhood data with, patch ID, coordinates (x,y) and species densities}

\item{disp_inds}{number immigrants to distribute (for each species)}

\item{delta}{coefficients determining strength of frequency-dependence}

\item{discrete}{whether or not immigrants should be distributed as discrete individuals}
}
\value{
a dataframe with neighborhood coordinates and immigrants densities (NOT updated densities!)
}
\description{
Applies weights to each neighborhood cell and distributes emigrants to each accordingly
}
