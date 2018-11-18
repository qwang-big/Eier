\name{edgeRank}
\alias{edgeRank}
\title{
Assign weights to edges
}
\description{
Assign weights to edges by summing the ranks of the two connected nodes
}
\usage{
edgeRank(pg, g, cutoff, min, reverse)
}

\arguments{

\item{pg}{
Vector contains sorted gene names. Highly ranked genes are ordered firstly if they are returned from \code{\link{pageRank}}.
}
\item{g}{
An igraph network object which consists the gene names. 
}
\item{min}{
Minimum scores assigned to the low ranking genes.
}
\item{cutoff}{
A numeric value between 0 and 1 indicating the percentage of genes are considered as low ranking / discarded.
}
\item{reverse}{
Logical indicating the order of the list. TRUE means highly ranked genes are ordered firstly, otherwise FALSE.
}
}

\details{
This function prepares the weighted network according to the gene ranks for further clustering.
}

\value{
An igraph network object with weights assigned to the edges.
}

\author{
Qi Wang
}
\examples{
data(hprd)
data(sample)
g <- edgeRank(res$pg, hprd, cutoff=0.9, min=1e-7)
}
