\name{exportJSONpathways}
\alias{exportJSONpathways}
\title{
Export Pathways
}
\description{
Export enriched pathways to a JSON file
}
\usage{
exportJSONpathways(enriched, prefix, n)
}

\arguments{

\item{enriched}{
Lists of data.frames exported by \code{\link{annotNets}}
}
\item{file}{
The prefix of the JSON filename to be exported.
}
\item{n}{
Number of pathways to be exported from each database, the pathways are sorted from smallest to the largest by adjusted p-values.
}
}

\details{
This functions writes a JSON file containing the pathways, p-values, and gene names in each pathway which can be further display with the web application \code{net}.  
}

\author{
Qi Wang
}
\examples{
data(sample)
exportJSONpathways(res$ga, 'CLL', n=15)
}
