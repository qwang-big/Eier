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
testBed <- tempfile()
control <- tempfile()
enhancer <- data.frame(chr="chr1", start=c(1,4,7,10), end=c(2,5,8,11), id=c("GH123456789","GH987654321","prom_1","prom_2"))
write.table(data.frame(chr="chr1", start=c(1,4,7), end=c(2,5,8)), testBed, col.names=F, row.names=F, quote=F, sep="\t")
write.table(data.frame(chr="chr1", start=c(1,7,10), end=c(2,8,11)), control, col.names=F, row.names=F, quote=F, sep="\t")
i <- filterPeak(c(testBed,control),enhancer,group=1:2,filter="enhancer")
}
