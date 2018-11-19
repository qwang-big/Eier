\name{writeMarkers}
\alias{writeMarkers}
\title{
Write Marker Genes to File
}
\description{
Write marker genes specific to each test case to file
}
\usage{
writeMarkers(id, id2, file)
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
genes=list(marker=c("BRAF","KRAS","MYC","TP53"))
writeMarkers(genes,"CLL")
}
