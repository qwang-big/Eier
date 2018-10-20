\name{exportJSONnets}
\alias{exportJSONnets}
\title{
Export Network
}
\description{
Export network to a JSON file
}
\usage{
exportJSONnets(gs, file)
}

\arguments{

\item{gs}{
Lists of igraph networks to export.
}
\item{file}{
The filename of the JSON file to be exported.
}
}

\details{
This functions writes a JSON file containing the highly-ranked subnetworks which can be further display with the web application \code{net}. 
}

\author{
Qi Wang
}
\examples{
exportJSONnets(gs, 'CLLnets.json')
}
