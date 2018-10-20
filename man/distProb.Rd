\name{distProb}
\alias{distProb}
\title{
Get Interaction Probability
}
\description{
Get interaction probability by power-transformation of promoter-enhancer distances
}
\usage{
distProb(file, model)
}

\arguments{

\item{file}{
Tab-delimited file contains three fields per line: (1) enhancer ID, (2) promoter ID, and (3) their distance (in base pairs, negative values indicate enhancer is in the upstream of the promoter). 
}

\item{model}{
Linear model created using \code{\link{genModel}}. 
}
}

\details{
The function takes in the coefficient (a) and intercept (b) from genModel, and applied to the presumed probability density function y=e^ax + b, to get the estimated probability for each given distances between the enhancer to the promoter. 
}

\value{
A data.frame with the last column converted the distances to relative probability.
}

\author{
Qi Wang
}

\seealso{
\code{\link{genModel}}
}
\examples{
}
