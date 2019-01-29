\name{readBW}
\alias{readBW}
\title{
Read data from one BigWig file
}
\description{
Compute sum score of BigWig file over an input genomic region
}
\usage{
readBW(filename, seqname, start, end)
}

\arguments{

\item{filename}{
filename of the BigWig file.
}
\item{seqname}{
chromosome ID.
}
\item{start}{
numeric value of start position.
}
\item{end}{
numeric value of end position.
}
}

\details{
'readBW' returns a numeric value which sums scores from the corresponding genomic loci of the BigWig file. 
}

\value{
A numeric value of the density scores of the BigWig file.
}

\author{
Qi Wang
}
\examples{
data <- readBW(system.file("data", "test.bw", package="irene"),"chr1",100000,123000)
}
