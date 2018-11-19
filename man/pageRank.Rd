\name{pageRank}
\alias{pageRank}
\title{
Rank meta-genes with the PageRank algorithm
}
\description{
Rank meta-genes according to the scores from the promoter and targeting enhancers
}
\usage{
pageRank(pcs, x, damping = 0.85, dWeight = 1e-99, fun = NULL,
    maxTry = 10, gene.id = TRUE, rewire = FALSE, statLog = NULL)
}

\arguments{

\item{pcs}{
data frame containing the ID of genomic regions, and a series of PCs for ranking. 
}
\item{pe}{
data frame representing the likelihood of enhancer-promoter interactions. The first column is the enhancer Id, and the second column is the promoter Id, and the third is the likelihood of enhancer-promoter interaction. 
}
\item{damping}{
the damping factor in the PageRank algorithm. 
}
\item{gene.id}{
logical, only output gene names if true, otherwise both gene names and their PageRank scores will be output. 
}
}

\details{
This function ranks meta-genes based on the 'personalized' PageRank algorithm in the 'igraph' implementation. The node weights are linearly transformed from the computed PCs, and edge weights are from the likelihood of enhancer-promoter interactions. The enhancer-promoter network are built from a bipartite directed graph, wherein the directions are from the enhancer nodes to their target promoter nodes. 
}

\value{
A sorted vector representing the PageRank order of meta-gene names, which highly ranked genes will appear on top. 
}

\author{
Qi Wang
}
\examples{
data(CLL)
data(sample)
j <- c(1,2,6,8)
res <- dPCA(meta, bed, data, datasets=j, transform=j, normlen=j, processedData=TRUE, verbose=TRUE)
pg <- pageRank(res$gr, h1)
}
