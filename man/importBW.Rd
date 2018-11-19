\name{importBW}
\alias{importBW}
\title{
Import density data from BigWigs
}
\description{
Compute sum score of BigWig files over input genomic regions
}
\usage{
importBW(meta, granges, threads)
}

\arguments{

\item{meta}{
a character vector of file names or a meta object containing a 'file' field.
}
\item{granges}{
genomic ranges in either data.frame or granges object. It is supposed to contain four columns, which represent chromosome name, start, end, and ID of the genomic region. 
}
\item{threads}{
number of cores to use.
}
}

\details{
'importBW' returns a numeric matrix with each row as the name of genomic loci, and each column as the BigWig files. The matrix contains sum scores of each BigWig file from corresponding genomic locus. 
}

\value{
A data frame containing a representation of the density scores of the BigWig files.
}

\author{
Qi Wang
}
\examples{
meta <- data.frame(file=c(system.file("data", "test.bw", package="irene"), 
	system.file("data", "ctrl.bw", package="irene")),group=1:2,dataset=1)
bed  <- data.frame(chr=c("chr1","chr1"),start=c(100000,187000),end=c(123000,188000),id=1:2)
data <- importBW(meta, bed)
}
