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
Selected database IDs, see \code{listEnrichrDbs} from \code{enrichR} for more information. By default, it uses KEGG (ID=94) and WikiPathways (ID=95) for enrichment. 
}
}

\value{
A list of data.frame containing enrichR outputs for each network
}

\author{
Qi Wang
}
\examples{
data(sample)
annotNets(cll$gs, sel=c(94,95))
}
