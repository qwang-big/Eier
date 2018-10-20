\name{writeData}
\alias{writeData}
\title{
Write Epigenetic Data
}
\description{
Write epigenetic intesity data to a file.
}
\usage{
writeData(data, groupLabels, filename, n, name, intTemp)
}

\arguments{

\item{data}{
A data.frame of intensity data to export. 
}
\item{groupLabels}{
Labels of biological data, use "Treament" or "Control" to indicate the group of data.
}
\item{filename}{
The filename to be saved as. 
}
\item{n}{
Integer of n-tiles for seperating the data. 
}
\item{name}{
Labels of experiment data, use "WGBS", "fractional" or "Bisulfite" to indicate DNA methylation data, H3KXXXX to indicate histone modification data.
}
\item{intTemp}{
Logical whether to write template files for holding the chromatin interation data. 
}
}

\details{
This function writes the epigenetic intensity data per loci sorted with chromesomal positions to a csv file. The data must include all the columns of the returned genomic regions object from \code{\link{dPCA}}
}


\author{
Qi Wang
}

\seealso{
\code{\link{dPCA}}
}

\examples{
writeData(cll$gr, c("CLL", "mCLL"), "CLL", intTemp=FALSE)
}
