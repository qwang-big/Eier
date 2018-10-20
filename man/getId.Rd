\name{getId}
\alias{getId}
\title{
Get Desired IDs
}
\description{
Get gene, enhancer or promoter IDs
}
\usage{
getId(id ,type)
}

\arguments{

\item{id}{
Vector of IDs to extract.
}
\item{type}{
Type of IDs to extract, must be one of the following values:
\itemize{
\item{promoter}{
Promoter ends with _1 or _2, ..., indicating the isotypes of promoter. 
}
\item{enhancer}{
Enhancer starts with "GH" (Genehancer IDs) or simply likes "chr1_999_1000".
}
\item{gene}{
Promoter IDs striped of suffix. 
}
\item{all}{
Both gene IDs and enhancer IDs.
}
}
}
}

\value{
Vector of desired IDs
}

\author{
Qi Wang
}
\examples{
promId <- getId(id, type="promoter")
}
