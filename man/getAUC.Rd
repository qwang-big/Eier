\name{getAUC}
\alias{getAUC}
\title{
Get AUC of ROC for Gene Ranks
}
\description{
Get AUC of ROC for gene ranks according the given list of marker genes
}
\usage{
getAUC(x, pos)
}

\arguments{

\item{x}{
A list of genes names sorted by the ranks.
}
\item{pos}{
Vector of the marker genes.
}
}

\details{
The function returns multiple AUCs according to the order of marker genes in the ranking list. 
}

\author{
Qi Wang
}
\examples{
data(markers)
data(sample)
getAUC(res$pg, markers$CLL)
}
