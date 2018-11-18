\name{exportJSONnets}
\alias{exportJSONnets}
\title{
Export Network
}
\description{
Export network to a JSON file
}
\usage{
exportJSONnets(gs, prefix)
}

\arguments{

\item{gs}{
Lists of igraph networks to export.
}
\item{file}{
The prefix of the JSON filename to be exported.
}
}

\details{
This functions writes a JSON file containing the highly-ranked subnetworks which can be further display with the web application \code{net}. 
}

\author{
Qi Wang
}
\examples{
data(sample)
exportJSONnets(res$gs, 'CLL')
}
