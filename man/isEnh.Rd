\name{isEnh}
\alias{isEnh}
\title{
Test an input ID
}
\description{
Test if the input ID is an enhancer ID
}
\usage{
isEnh(id)
}

\arguments{

\item{id}{
the ID to be tested
}
}

\details{
Enhancer ID starts with "GH" (Genehancer symbol and is more than ten characters) or simply likes "chr7_74540885_74541801".
}

\value{
logical indicating whether the input ID is an enhancer ID
}

\author{
Qi Wang
}
\examples{
isEnh("GH08I085513")
isEnh("GHR_1")
isEnh("chr7_74540885_74541801")
}
