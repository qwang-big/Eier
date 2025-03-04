\name{exportD}
\alias{exportD}
\title{
Export Dobs
}
\description{
Export Dobs to a csv file
}
\usage{
exportD(Dobs, labels, prefix)
}

\arguments{

\item{Dobs}{
The D matrix, which contains the observed differences between the two conditions returned by \code{\link{dPCA}}.
}
\item{labels}{
Vector consists of the labels of each component in the D matrix. 
}
\item{file}{
The filename of the csv file to be exported.
}
}

\details{
This functions writes a csv file containing the composition of the components which can be further display with the web application \code{net}, \code{rank}, and \code{rna}.  
}

\author{
Qi Wang
}
\examples{
data(sample)
exportD(res$Dobs, c('K4me1','K4me3','K27ac','Meth'), "CLLcomp.csv")
}
