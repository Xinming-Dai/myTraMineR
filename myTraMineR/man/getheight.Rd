% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getheight.R
\name{getheight}
\alias{getheight}
\title{Get the cut height for cutting a hierarchical clustering tree}
\usage{
getheight(seqdistance, cmethod = "complete", idxs = 1:10, fq)
}
\arguments{
\item{seqdistance}{A distance matrix or a distance array returned by TraMineR::seqdist() function.}

\item{cmethod}{This argument will be passed to \code{hclust()}. It takes the same values as \code{method} argument in \code{hclust()}. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC). The default clustering method is "complete", complete-linkage clustering.}

\item{idxs}{A integer or an array of integers. The Default is 1:10, meaning replacing sequences in the 10 largest clusters with their representatives. If idxs=0, then all clusters will be replaces.}

\item{fq}{A numeric number within the range [0,1].}
}
\value{
The cut height—a numeric number and the frequency, meaning how much these clusters will represent for the total.
}
\description{
This function returns a cut height, denoted as \code{h}, on a dendrogram using a specified clustering method (\code{cmethod}). If the dataset is clustered based on \code{h}, the first few most frequent groups, specified by \code{idxs}, will yield at least \code{fq} of the total sequences.
}
\examples{
getheight(seqdistance, fq=0.3)
getheight(seqdistance, fq=0.3, idxs=0)
}
