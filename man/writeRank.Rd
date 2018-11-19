\name{writeRank}
\alias{writeRank}
\title{
Write Gene Ranks to File
}
\description{
Write gene ranks derived from enhancers or promoters to file
}
\usage{
writeRank(id, id2, prefix, header = c("PromEnh", "PromOnly"))
}

\arguments{

\item{id}{
The sorted IDs of gene ranks derived both enhancers and promoters.
}
\item{id2}{
The sorted IDs of gene ranks derived promoters only.
}
\item{file}{
The filename to save as.
}
}

\details{
The function save two columns of gene ranks pair-wiselly to a csv file. 
}

\author{
Qi Wang
}
\examples{
data(sample)
writeRank(res$pg[[1]], res$pg$prom, "CLLrank.csv")
}
