\name{dPCA}
\alias{dPCA}
\title{
Differential Principal Component Analysis
}
\description{
Run differential principal component analysis with incorporated dPCA program
}
\usage{
dPCA(meta, bed, data, sampleId = NULL, groups = 1:2,
    datasets = NULL, transform = NULL, normlen = NULL,
    minlen = 50, lambda = 0.15, fun = function(x)
    sqrt(mean(x^2)), datasetLabels = NULL, groupLabels =
    NULL, qnormalize = TRUE, qnormalizeFirst = FALSE,
    normalize = FALSE, verbose = FALSE, interactive =
    FALSE, useSVD = FALSE, saveFile = FALSE, processedData
    = TRUE, removeLowCoverageChIPseq = FALSE,
    removeLowCoverageChIPseqProbs = 0.1, dPCsigns = NULL,
    nPaired = 0, nTransform = 0, nColMeanCent = 0,
    nColStand = 0, nMColMeanCent = 1, nMColStand = 0,
    dSNRCut = 5, nUsedPCAZ = 0, nUseRB = 0, dPeakFDRCut =
    0.5)
}

\arguments{

\item{meta}{
Meta information consists three columns (1) filenames, (2) biological groups, and (3) dataset IDs. Can be either data.frame or matrix. 
}
\item{bed}{
Input genomic regions in data.frame or matrix format, which needs to have four columns in exact this order: chromosome names, start position, end position, ID. 
}
\item{data}{
Data.frame or matrix which contains the intensity data for each input genomic regions. Columns are sorted in the order of filenames of meta information, can be produced using \code{\link{importBW}}. 
}
\item{sampleId}{
Vector of sample IDs used to be tested, equivalent to the line numbers of the meta information. 
}
\item{groups}{
Vector or List of group IDs used to be tested, specified in the group fields of the meta information. If the input is a list, for example \code{list(1:3,4)}, all IDs in the first vector from the list are considered as the same group. 
}
\item{datasets}{
Vector of dataset IDs used to be tested, specified in the datasets fields of the meta information. 
}
\item{transform}{
Vector of dataset IDs which need to be transformed, power-transformations are applied according to the lambda estimated by \code{boxcox} function. If lambda is equal to or less than 0, log-transformations are applied instead. 
}
\item{normlen}{
Vector of datasets IDs which need to normalized according to the length of the genomic regions. 
}
\item{minlen}{
Numeric value of minimum length of the genomic regions to be tested, genomic regions less than this threshold are discarded.
}
\item{normalize}{
Logical whether data will be normalized to the total library size, default is FALSE. 
}
\item{qnormalize}{
Logical whether quantile will be applied to all samples, default is TRUE. 
}
\item{fun}{
Transformation function to be applied to PCs, which produces another PCx field in the output incorporating one or several PCs. Here are some examples:
\itemize{
\item function(x) sqrt(mean(x^2))
\item function(x) log2(mean(2^x))
\item function(x) x[1]
\item function(x) abs(x[1])
\item function(x) x[1]-x[2]
\item function(x) x[1]-x[2]-x[3]
\item function(x) x[1]+x[2]
\item function(x) sum(x)
\item function(x) sum(abs(x))
\item function(x) x[1]/x[2]
\item function(x) max(abs(x))
\item function(x) x[which.max(abs(x))]
}
and by default, function(x) sqrt(mean(x^2)) is used.
}
\item{verbose}{
Logical whether additional information and figures are shown during power-transformation and normalization. 
}
\item{processedData}{
Logical whether the processed data will be returned with other outputs. 
}
}

\details{
This function filters, normalizes and transforms the desired groups and datasets of the data, then forwards the processed data to a incorporated C program called dPCA (see PMID: 23569280), and outputs the dPCs and other information from the program. 
}

\value{
A named list contains three data.frames:
\itemize{
\item{gr}{
The predefined regions followed by the computed dPCs from dPCA, as well as PCx produced by the transformation function. If \code{processedData} is TRUE, the processed data will also be append to the data.frame. 
}
\item{Dobs}{
The D matrix, which contains the observed differences between the two conditions. This is the data analyzed by dPCA.
}
\item{proj}{
Estimated beta coefficients for each dPC.
}
}
}

\author{
Qi Wang
}
\examples{
data(CLL)
j <- c(1,2,6,8)
res <- dPCA(meta, bed, data, groups=1:2, datasets=j, transform=j, normlen=j, processedData=TRUE, verbose=TRUE)
}
