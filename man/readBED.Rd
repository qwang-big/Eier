\name{readBED}
\alias{readBED}
\title{
Read a BED File
}
\description{
Read a BED File into a data.frame
}
\usage{
readBED(file)
}

\arguments{

\item{file}{
The filename of BED file to read. 
}
}

\details{
The function reads a BED file, and assign the first column to chromosome name, second column to chromosomal start position, third column to chromosomal end position, and fourth column (if available) to loci ID.
}

\value{
Data.frame with chromosome name, start position, end position, and loci ID.
}

\author{
Qi Wang
}
\examples{
promoter <- readBED('Promoter.bed')
}
