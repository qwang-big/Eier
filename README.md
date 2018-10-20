# Eier: Enhancer-integrated epigenetic ranking for cancer and development

![logo](https://raw.githubusercontent.com/qwang-big/eier/master/images/Eier.png)

[![Build Status](https://travis-ci.org/qwang-big/Eier.svg?branch=master)](https://travis-ci.org/qwang-big/Eier)

*Eier* is an R package which allows you
  - Find significantly altered genes between two biological conditions from histone ChIP-Seq and DNA methylation tracks (BigWig)
  - Focus on the genes which are associated with more extensive epigenetic modification over targetting enhancers.
  - Annotate and visulise the global epigenetic variances with network anlaysis.

# INTRODUCTION
*Eier* is developed for two purposes in epigenetic ranking:
  - Integrate several epigenetic marks
  - Incorporate enhancers

The whole idea has been demostrated in the Chapter III of Qi Wang's disseration. For the above purposes, we employed singular value decomposition [dPCA](http://www.biostat.jhsph.edu/~hji/dpca/) and random walk ranking [PageRank](http://igraph.org/r/doc/page_rank.html). The epigenetic alterations over genomic regulatory elements between two groups are presented as dPC scores and further to PageRank scores during solving the enhancer-promoter relationships. To begin with, *Eier* starts from the following data inputs. 

# Prerequisites

## Genomic regions to be tested
User need to provide promoters and enhancers regions so that comparisons for differential epigenetic modification can be made for those regions. We also provide pre-defined promoters from [The Eukaryotic Promoter Database](http://epd.vital-it.ch/index.php) and pre-defined from [GeneHancer] (http://www.genecards.org/Guide/GeneCard). 
Given that the enhancers from GeneHancer database are an ensemble of all tissues/cell types, one may need to filter out the unspecific ones by overlapping with the enhancer-specific histone marks, e.g., histone H3 lysine 4 monomethylation (H3K4me1) or the H3K27 acetylation (H3K27ac) marks. 

## Promoter-enhancer (P-E) interactions
Enhancer within 1Mb distance to the TSS are considered potential interating region. In a common sense, the interactions should not cross the topologically associating domains (TAD) boundary. Therefore, we compiled a cell-type specific enhancer-promoter interaction list by excluding the interactions not within the same TAD from [GSE87112](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87112). Enhancer-promoter interactions probability are estimated using a power-law decay function based on the distances to the TSS. Enhancer-promoter distances for sample tissues can be downloaded from [https://github.com/qwang-big/crl-data/PEdistances](https://github.com/qwang-big/crl-data/tree/master/PEdistances), including:
 - H1: H1 e

In practice, as TADs between different cell types are relative conserved ([Schmitt AD, 2016](#schmitt-ad)), one can use the H1 cell line in case the TAD of corresponding cell type is not available. 

First Tab:
```sh
$ node app
```


## Epigenetic intensity data
*Eier* requires user to provide [**BigWig**](https://genome.ucsc.edu/goldenpath/help/bigWig.html) format to represent the sequencing density of BS-Seq or ChIP-Seq. User need to create a **data.frame** to indicate the location of the BigWig files, as well as groups and experiment types (as *dataset* column). Here is a sample as follows:

|    | file                                                                               | group | dataset |
|----|------------------------------------------------------------------------------------|-------|---------|
| 1  | UCSD.H1_Derived_Mesenchymal_Stem_Cells.Bisulfite-Seq.methylC-seq_h1-msc_r1a.wig.bw | 1     | 1       |
| 2  | UCSD.H1_Derived_Mesenchymal_Stem_Cells.Bisulfite-Seq.methylC-seq_h1-msc_r2a.wig.bw | 1     | 1       |
| 3  | UCSD.H1_Derived_Mesenchymal_Stem_Cells.H3K27ac.SK436.wig.bw                        | 1     | 2       |
| 4  | UCSD.H1_Derived_Mesenchymal_Stem_Cells.H3K27ac.SK438.wig.bw                        | 1     | 2       |
| 5  | UCSD.H1_Derived_Mesenchymal_Stem_Cells.H3K27me3.SK437.wig.bw                       | 1     | 3       |
| 6  | UCSD.H1_Derived_Mesenchymal_Stem_Cells.H3K27me3.SK439.wig.bw                       | 1     | 3       |
| 7  | UCSD.H1_Derived_Mesenchymal_Stem_Cells.H3K9ac.SK518.wig.bw                         | 1     | 4       |
| 8  | UCSD.H1_Derived_Mesenchymal_Stem_Cells.H3K9ac.SK519.wig.bw                         | 1     | 4       |
| 9  | UCSD.H1_Derived_Mesenchymal_Stem_Cells.H3K9me3.SK507.wig.bw                        | 1     | 5       |
| 10 | UCSD.H1_Derived_Mesenchymal_Stem_Cells.H3K9me3.SK508.wig.bw                        | 1     | 5       |
| 11 | UCSD.H1_Derived_Mesenchymal_Stem_Cells.Input.SK443.wig.bw                          | 1     | 6       |
| 12 | UCSD.H1_Derived_Mesenchymal_Stem_Cells.Input.SK444.wig.bw                          | 1     | 6       |
| 13 | GSM429321_UCSD.H1.Bisulfite-Seq.methylC-seq_h1_r2a.bw                              | 2     | 1       |
| 14 | GSM429322_UCSD.H1.Bisulfite-Seq.methylC-seq_h1_r2b.bw                              | 2     | 1       |
| 15 | GSM432685_UCSD.H1.Bisulfite-Seq.methylC-seq_h1_r1a.bw                              | 2     | 1       |
| 16 | GSM432686_UCSD.H1.Bisulfite-Seq.methylC-seq_h1_r1b.bw                              | 2     | 1       |
| 17 | UCSD.H1.H3K27ac.LL313.bw                                                           | 2     | 2       |
| 18 | UCSD.H1.H3K27ac.SAK270.bw                                                          | 2     | 2       |
| 19 | UCSD.H1.H3K27me3.LL241.bw                                                          | 2     | 3       |
| 20 | UCSD.H1.H3K27me3.LL314.bw                                                          | 2     | 3       |
| 21 | UCSD.H1.H3K27me3.YL95.bw                                                           | 2     | 3       |
| 22 | UCSD.H1.H3K9ac.LL240.wig.bw                                                        | 2     | 4       |
| 23 | UCSD.H1.H3K9ac.SAK68.wig.bw                                                        | 2     | 4       |
| 24 | UCSD.H1.H3K9me3.AK54.wig.bw                                                        | 2     | 5       |
| 25 | UCSD.H1.H3K9me3.LL218.wig.bw                                                       | 2     | 5       |
| 26 | UCSD.H1.H3K9me3.YL75.wig.bw                                                        | 2     | 5       |
| 27 | UCSD.H1.H3K9me3.YL77.wig.bw                                                        | 2     | 5       |
| 28 | UCSD.H1.Input.AK57.wig.bw                                                          | 2     | 6       |
| 29 | UCSD.H1.Input.DM219.wig.bw                                                         | 2     | 6       |
| 30 | UCSD.H1.Input.LL-H1-I1.wig.bw                                                      | 2     | 6       |
| 31 | UCSD.H1.Input.LL-H1-I2.wig.bw                                                      | 2     | 6       |
| 32 | UCSD.H1.Input.LLH1U.wig.bw                                                         | 2     | 6       |
| 33 | UCSD.H1.Input.YL154.wig.bw                                                         | 2     | 6       |
| 34 | UCSD.H1.Input.YL208.wig.bw                                                         | 2     | 6       |
| 35 | UCSD.H1.Input.YL262.wig.bw                                                         | 2     | 6       |
| 36 | UCSD.H1.Input.YL328.wig.bw                                                         | 2     | 6       |

The **BigWig** files for our test cases were converted from wig files retrieved from [NIH Roadmap Epigenomics Project data gateway](https://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/), BLUEPRINT, and [CEEHRC](http://www.epigenomes.ca/), respectively.  Precompiled R objects in our test cases can be retrieved from  

# Usage
## R environment for the test cases.
The following options were set up during performing our tests. 
```r
options(stringsAsFactors = FALSE)
```

## Load pre-defined regions
The pre-defined promoters and enhancers regions with corresponding IDs (precompiled [hg19](https://github.com/qwang-big/crl-data/blob/master/promenh.hg19.bed), [hg38](https://github.com/qwang-big/crl-data/blob/master/promenh.hg38.bed)) need to be loaded first with: 
```r
bed <- read.table('promenh.hg19.bed')
```

## Import sequencing density data from BigWig files
Read in the table containing file location, group, dataset information as a **data.frame** (named *meta* in the following case), the use 
```r
data <- importBW(meta, bed)
```
to import the density for each regions. *ImportBW* use an external C function from [libBigWig](), which was compiled with GNU gcc and linked with a static C library during *Eier* installation on Linux.  For Windows/Mac OS users, you need to create your platform specific static library (**libBigWig.a**, under *src/libBigWig*) by your own. 

## Filter out unspecific enhancers
We only took the regions which are more likely to be true enhancers, therefore we use the following function to get the indices which overlapped with enhancer histone marks (H3K4me1 in the following case). The peaks identified by Roadmap Epigenetics Project were retrieved from [http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated](http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated), and run:
```r
i <- filterPeak(c("E003-H3K4me1.narrowPeak","E006-H3K4me1.narrowPeak"), bed, group=c(1,2))
```
**Precompiled R data files in our test cases () contain all the above mentioned objects (meta, bed, data, i) for the following tests. **

## Measure combinatorial effect of epigenetic alterations
We use [dPCA](http://www.biostat.jhsph.edu/~hji/dpca/) to measure combinatorial effect of epigenetic alterations. The software is already integrated into *Eier* as an external C function, which was compiled with GNU gcc and linked with a static C library during *Eier* installation on Linux.  For Windows/Mac OS users, you need to create your platform specific static library (**libdpca.a**, under *src/dpca*) by your own. 
```r
res <- dPCA(meta, bed[i,], data[i,], datasets=1:3, transform=1:3, normlen=1:3, verbose=T)
```
The output *res* is a named list which has three keys: *gr*, *Dobs*, and *proj*.

* **gr**
 contains the pre-defined regions followed by the computed PCs from dPCA. A sample output may look like this:

| seqnames | start     | end       | id          | PC1      | PC2      | PC3      |
|---------------|------------|------------|-------|-------------|----------|----------|
|     chr5 | 159099769 | 159099829   |   EBF1_1 | 22.12935 | -11.913720 | -1.9786781 |
|   chr19 | 5386940 1 |  53869461   |  MYADM_3 | 22.12694 | -8.876568 | 3.7468742  |
|     chr4 | 52862267 | 52862327  |  RASL11B_1 | 21.89510 | -14.516893 | -4.8862062 |
|     chr8 | 141099562 | 141099711 | GH08I141099 | 21.88550  |  7.578159 | 1.3920079 |
|     chr6 | 131063357 | 131063417  | EPB41L2_4 | 21.57229 | -12.886749 | -0.8408109 |
|     chr2 | 158457054 | 158457114   |   PKP4_1 | 21.42755 | -8.421447 | -4.6387676 |

The table is sorted descendingly by PC1

* **Dobs**
The D matrix, which contains the observed differences between the two conditions. This is the data analyzed by dPCA.

* **proj**
Estimated beta coefficients for each dPC.

## Infer promoter-enhancer interaction probabilities according to distances
Use the following function to convert the promoter-enhancer interactions in a given range of TSS (1Mb in the data provided) to probabilities of interaction.

* Load P-E interactions in H1ESC cell lines, the three columns in *H1* data.frame are enhancer IDs, promoter IDs, P-E distances, respectively. 
```r
H1 <- read.table('H1.hg19.pair')
```

* Transform the P-E interactions from bp to Mb:
```r
H1[,3]=abs(H1[,3]/1000000)
```

* Apply a power-decay function to represent the likelihood of P-E interactions. The choice of power-coefficient is explained in Figure 3.14b of Qi Wang's dissertation. We recommand using -20 for all the test cases. 
```r
H1[,3]=exp(-20*H1[,3]+1)
```

## Get rank of the promoters, which explain the significance of promoter alterations
Use the following function to get the order of promoters sorted by _nPC_, the ones on top of the list indicates higher importance of the corresponding gene. 
```{r, eval = FALSE}
id <- getId(res$gr$id, type="gene")
```

## Using PageRank to rank the importance of the alterations.
The PageRank will take in the alteration contribution from enhancers, resulting the ranks of associated promoters get promoted if targeted by high-rank enhancers. 
```{r, eval = FALSE}
p <- pageRank(res$gr$id, P, fun="logit")
```

## Use literature derived marker genes as a metric to evaluate the ranking
```{r, eval = FALSE}
rid <- getId(names(p), type="gene")
MSC <- unique(read.table('~/tmp/pk/msc', stringsAsFactors = FALSE)[,1])
plotRank(list('promoter'=id, 'enhancer'=rid), MSC)
```


# Case studies
Epigenetic and expression data were downloaded from [**EdaccData Release-9**](http://genboree.org/EdaccData/Release-9/experiment-sample/), [**CEEHRC**](http://www.epigenomes.ca/), [**IHEC**](http://ihec-epigenomes.org/). 

# Acknowledgements
The results presented here are in part based upon data generated by The Canadian Epigenetics, Epigenomics, Environment and Health Research Consortium (CEEHRC) initiative funded by the Canadian Institutes of Health Research (CIHR), Genome BC, and Genome Quebec. Information about CEEHRC and the participating investigators and institutions can be found at http://www.cihr-irsc.gc.ca/e/43734.html. 

# References
(### schmitt-ad)
Schmitt AD, Hu M, Jung I, et al. A Compendium of Chromatin Contact Maps Reveals Spatially Active Regions in the Human Genome. Cell Reports. 2016;17(8):2042–2059.

License
----

MIT

# Session Info
```r
pander(sessionInfo(), compact = FALSE)

```

**R version 3.2.2 (2015-08-14)**

**Platform:** x86_64-pc-linux-gnu (64-bit) 

**locale:**
_LC_CTYPE=en_GB.UTF-8_, _LC_NUMERIC=C_, _LC_TIME=en_GB.UTF-8_, _LC_COLLATE=en_GB.UTF-8_, _LC_MONETARY=en_GB.UTF-8_, _LC_MESSAGES=en_GB.UTF-8_, _LC_PAPER=en_GB.UTF-8_, _LC_NAME=C_, _LC_ADDRESS=C_, _LC_TELEPHONE=C_, _LC_MEASUREMENT=en_GB.UTF-8_ and _LC_IDENTIFICATION=C_

**attached base packages:** 

* parallel 
* stats4 
* stats 
* graphics 
* grDevices 
* utils 
* datasets 
* methods 
* base 


**other attached packages:** 

* pander(v.0.6.0) 
* BiSeq(v.1.11.0) 
* Formula(v.1.2-1) 
* SummarizedExperiment(v.1.0.1) 
* Biobase(v.2.30.0) 
* GenomicRanges(v.1.22.4) 
* GenomeInfoDb(v.1.6.1) 
* IRanges(v.2.4.8) 
* S4Vectors(v.0.8.11) 
* BiocGenerics(v.0.16.1) ​

