\name{combineGrpRgn}
\alias{combineGrpRgn}
\title{
Intersect Combined Regions
}
\description{
Combined genomic regions for each group, then get the intersection
}
\usage{
combineGrpRgn(files, grp)
}

\arguments{

\item{files}{
Genomic region filenames in BED format.
}
\item{grp}{
Groups vector in correspondance to the files.
}
}

\details{
This function calls combineRgn for each group internally, and intersect the combined genomic regions.
}

\value{
Genomic ranges of the combined regions.
}

\author{
Qi Wang
}
\examples{
}
