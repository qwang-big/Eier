\name{plotD}
\alias{plotD}
\title{
Plot D Matrix
}
\description{
Plot the principal components of D matrix analyzed by dPCA.
}
\usage{
plotD(Dobs, labels, captions, title, saveFile)
}

\arguments{

\item{data}{
D matrix to be analyzed.
}
\item{labels}{
Vector of the labels of experiment.
}
\item{captions}{
Vector of the titles of principal components.
}
\item{title}{
Title appeared on the plot.
}
\item{saveFile}{
Logical whether to export the plot as a png file
}
}

\details{
This function performs a principal components analysis on the D matrix, and plot the components per sample. 
}

\seealso{
\code{\link{dPCA}}
}

\author{
Qi Wang
}
\examples{
data(sample)
plotD(res$Dobs, c('K4me1','K4me3','K27ac','Meth'), scales::percent(res$proj[2,]))
}
