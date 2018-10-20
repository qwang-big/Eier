\name{filterPeak}
\alias{filterPeak}
\title{
Filter Peak Sets
}
\description{
Filter peak sets with reference genomic regions
}
\usage{
filterPeak(files, ref, group, filter)
}

\arguments{

\item{files}{
Vector of filenames of the peak sets used to filter the reference genomic regions.
}
\item{ref}{
Data.frame of reference genomic regions needed to be filtered.
}
\item{group}{
Vector of groups used to combine the regions, this must be of the same length to the vector of filenames. If groups is set to NULL, \code{\link{combineRgn}} will be called instead of \code{\link{combineGrpRgn}}. 
}
\item{filter}{
Vector of which type of genomic regions to filter, the value must be one of "enhancer", "promoter" or "both".
}
}

\details{
This function combines peak sets with or without group information, then \code{findOverlaps} will be called to get IDs of the overlaps, and intersect with the IDs of the reference genomic regions. 
}

\value{
Vector of indices of overlapped reference genomic regions. 
}

\author{
Qi Wang
}
\seealso{
\code{\link{combineRgn}}, \code{\link{combineGrpRgn}}, \code{\link{filterRgn}}
}
\examples{
enhancer <- filterPeak(c("mCLL.bed","CLL.bed"), enhancer)
}
