\name{plotRank}
\alias{plotRank}
\title{
Plot ROC for Gene Ranks
}
\description{
Plot ROC for gene ranks according the given list of marker genes
}
\usage{
plotRank(x, pos, title, file, perc, captions)
}

\arguments{

\item{x}{
A list of genes names sorted by the ranks.
}
\item{pos}{
Vector of the marker genes.
}
\item{title}{
Title of the plot.
}
\item{file}{
The filename to be saved as. 
}
\item{perc}{
Numeric value between 0 and 1 indicating how much fraction of marker genes to use. 
}
\item{captions}{
Vector showing legend captions of the plot. 
}
}

\details{
The function plots multiple recovery curves according to the order of marker genes in the ranking list. 
}

\author{
Qi Wang
}
\examples{
data(markers)
data(sample)
plotRank(res$pg, markers$CLL)
}
