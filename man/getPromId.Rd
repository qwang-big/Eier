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
data(sample)
res$pg$prom <- getPromId(res$gr)
}
