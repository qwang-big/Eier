\name{genModel}
\alias{genModel}
\title{
Generate Probability Model
}
\description{
Create probability density function from promoter-contacts counts
}
\usage{
genModel(file, minDist = 0, minCount = 3, title = "", verbose = TRUE) 
}

\arguments{

\item{file}{
A tab-delimited file. First column contains distances to the promoters, second column contains counts of the contacts.
}
\item{minDist}{
Minimum distances required to infer the probability density function.
}
\item{minCount}{
Minimum contact counts required to infer the probability density function.
}
\item{title}{
Title to be appeared in the probability density plot.
}
\item{verbose}{
Logical which controls whether to show the probability density plot and additional information.
}
}

\details{
The function tries to fit the counts to a presumed exponential function y=e^ax + b, which x is the distances to the promoters.
}

\value{
A fitted model object
}

\author{
Qi Wang
}
\examples{
model = genModel(file.path(system.file("data", package="crl"), "CD34.pd"), minCount=10)
}
\seealso{
\code{\link{distProb}}
}
