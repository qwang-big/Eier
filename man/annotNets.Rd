\name{annotNets}
\alias{annotNets}
\title{
Annotate Networks
}
\description{
Annotate lists of igraph networks to pathways using enrichR
}
\usage{
annotNets(gs, sel)
}

\arguments{

\item{gs}{
Lists of igraph networks to annotate.
}

\item{sel}{
Selected database IDs, see \code{listEnrichrDbs} from \code{enrichR} for more information.
}
}

\value{
A list of data.frame containing enrichR outputs for each network
}

\author{
Qi Wang
}
\examples{
annotNets(gs, sel=c(94,95))
}
