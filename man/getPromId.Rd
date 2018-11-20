\name{getPromId}
\alias{getPromId}
\title{
Get Promoter IDs
}
\description{
Get promoter IDs sorted by PC from a data.frame of genomic regions
}
\usage{
getPromId(gr ,pc)
}

\arguments{

\item{gr}{
data.frame of genomic regions output by \code{dPCA}, which contains the id and PC
}
\item{pc}{
the name of the PC used to sort the promoter IDs, must be contained in the column names of the data.frame
}
}

\value{
Vector of sorted promoter IDs
}

\author{
Qi Wang
}
\examples{
data(markers)
data(CLL)
j <- c(1,2,6,8)
res <- dPCA(meta, bed, data, datasets=j, transform=j, normlen=j, processedData=TRUE, verbose=TRUE)
prom <- getPromId(res$gr)
}
