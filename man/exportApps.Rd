\name{exportApps}
\alias{exportApps}
\title{
Export Apps
}
\description{
Export web applications for visualizing the results
}
\usage{
exportApps(name, markers, exdir = ".", check = TRUE)
}

\arguments{

\item{name}{
The prefix of the experiement name. All the output filenames are followed by the prefix. 
}
\item{genes}{
A list of gene names. The names of the list are used to present the source of marker genes. 
}
\item{exdir}{
The directory name for writing the files to. 
}
\item{check}{
logical for checking the dependencies.
}
}

\details{
This functions writes HTML files for visualizing the results from previous outputs.  
}

\author{
Qi Wang
}
\examples{
genes=list(marker=c("BRAF","KRAS","MYC","TP53"))
exportApps("CLL", genes, check=FALSE)
}
