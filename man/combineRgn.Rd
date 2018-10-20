\name{combineRgn}
\alias{combineRgn}
\title{
Combine Genomic Regions
}
\description{
Combine genomic regions from BED files
}
\usage{
combineRgn(files)
}

\arguments{

\item{files}{
Genomic region filenames in BED format.
}
}

\details{
This function computes the union of genomic regions if input files are more than one, otherwise will return the genomic regions in the only file.
}

\value{
Genomic ranges of the combined regions.
}

\author{
Qi Wang
}
\examples{
}
