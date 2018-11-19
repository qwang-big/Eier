\name{exportMultinets}
\alias{exportMultinets}
\title{
Export Network Communities
}
\description{
Export network communities according to clustering of edge weights
}
\usage{
exportMultinets(g, n = 5, steps = 4, rewire = FALSE, simplify = TRUE)
}

\arguments{

\item{g}{
The reference network with weights assigned to edges.
}
\item{n}{
Number of subnetworks to take in the initial partition.
}
\item{steps}{
The length of the random walks to perform.
}
\item{rewire}{
Logical whether rewiring edges of the graph is needed, set TRUE only if you want to validate the robustness of the communities.
}
}

\details{
The network must have weights assigned to the edges using \code{\link{edgeRank}} before running this function. The function selected top n ranking communities in the initial partition, and further breaks up the largest three network with another round of random walks. The subnetworks are sorted by the mean weights of all the nodes.
}

\value{
A list of subnetworks sorted by the mean weights of all the nodes.
}

\author{
Qi Wang
}

\seealso{
 \code{\link{edgeRank}}
}

\examples{
data(hprd)
data(sample)
g <- edgeRank(res$pg, hprd, cutoff=0.9, min=1e-7)
gs <- exportMultinets(g, n = 15, rewire = FALSE)
}
