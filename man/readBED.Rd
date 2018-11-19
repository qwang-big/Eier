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
#create a dummy bed file for testing
tf <- tempfile()
write.table(data.frame(chr=paste0("chr",1:5), start=1:5, end=6:10, id=paste0("gene",1:5)), tf, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
#read the bed file just created 
promoter <- readBED(tf)
}
