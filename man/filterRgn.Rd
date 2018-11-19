\name{filterRgn}
\alias{filterRgn}
\title{
Filter Genomic Regions
}
\description{
Filter genomic regions with reference genomic regions
}
\usage{
filterRgn(region, ref)
}

\arguments{

\item{region}{
Data.frame of input genomic regions.
}
\item{ref}{
Data.frame of reference genomic regions needed to be filtered.
}
}

\details{
This function calls \code{findOverlaps} and select the first entry from reference genomic regions if multiple input genomic regions overlapped with the reference genomic regions. 
}

\value{
Entries in the reference genomic regions overlapped with the input. 
}

\author{
Qi Wang
}
\seealso{
\code{\link{filterRgn}}
}
\examples{
testBed <- data.frame(chr="chr1", start=c(1,4,7), end=c(2,5,8))
control <- data.frame(chr="chr1", start=c(1,7,10), end=c(2,8,11))
filterRgn(makeGRangesFromDataFrame(testBed),makeGRangesFromDataFrame(control))
}
