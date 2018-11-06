# Irene: Integrative Ranking with Epigenetic Network of Enhancers

![logo](https://raw.githubusercontent.com/qwang-big/irene/master/images/Irene.png)

[![Build Status](https://travis-ci.org/qwang-big/irene.svg?branch=master)](https://travis-ci.org/qwang-big/irene)

*Irene* is an R package which allows you
  - Find significantly altered genes between two biological conditions from histone ChIP-Seq and DNA methylation tracks (BigWig)
  - Focus on the genes which are associated with more extensive epigenetic modification over targetting enhancers.
  - Annotate and visualize the global epigenetic variances with network analysis.

# INTRODUCTION
*Irene* is developed for two purposes in epigenetic ranking:
  - Integrate several epigenetic marks
  - Incorporate enhancers

The whole idea has been demonstrated in the Chapter III of Qi Wang's [dissertation](https://github.com/qwang-big/irene-web/blob/master/docs/chapter3.pdf). For the above purposes, we employed singular value decomposition [dPCA](http://www.biostat.jhsph.edu/~hji/dpca/) and random walk ranking [PageRank](http://igraph.org/r/doc/page_rank.html). The epigenetic alterations over genomic regulatory elements between two groups are presented as dPC scores and further to PageRank scores during solving the enhancer-promoter relationships. To begin with, *Irene* starts from the following data inputs. 

# Installation
*Irene* can be installed directly from GitHub with the help of *devtools* package:
```r
library(devtools)
install_github("qwang-big/irene")
```

# Prerequisites

## Genomic regions to be tested
User need to provide promoters and enhancers regions so that comparisons for differential epigenetic modification can be made for those regions. We also provide pre-defined promoters from [The Eukaryotic Promoter Database](http://epd.vital-it.ch/index.php) and pre-defined from [GeneHancer](http://www.genecards.org/Guide/GeneCard). 
Given that the enhancers from GeneHancer database are an ensemble of all tissues/cell types, one may need to filter out the unspecific ones by overlapping with the enhancer-specific histone marks, e.g., histone H3 lysine 4 monomethylation (H3K4me1) or the H3K27 acetylation (H3K27ac) marks. 

## Promoter-enhancer (P-E) interactions
Enhancer within 1Mb distance to the TSS are considered potential interacting region. In a common sense, the interactions should not cross the topologically associating domains (TAD) boundary. Therefore, we compiled a cell-type specific enhancer-promoter interaction list by excluding the interactions not within the same TAD from [GSE87112](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87112). Enhancer-promoter interactions probability are estimated using a power-law decay function based on the distances to the TSS. Enhancer-promoter distances for sample tissues can be downloaded from [https://github.com/qwang-big/irene-data/PEdistances](https://github.com/qwang-big/irene-data/tree/master/PEdistances), including:
 - [H1](https://github.com/qwang-big/irene-data/blob/master/PEdistances/H1.hg19.pair.gz): H1 human embryonic stem cell line
 - [MES](https://github.com/qwang-big/irene-data/blob/master/PEdistances/MES.hg19.pair.gz): H1 BMP4 derived mesendoderm cultured cells
 - [MSC](https://github.com/qwang-big/irene-data/blob/master/PEdistances/MSC.hg19.pair.gz): H1 derived mesenchymal stem cells
 - [NPC](https://github.com/qwang-big/irene-data/blob/master/PEdistances/NPC.hg19.pair.gz): H1 derived neural precursor cells
 - [TPC](https://github.com/qwang-big/irene-data/blob/master/PEdistances/TPC.hg19.pair.gz): H1 derived trophoblast stem cells

In practice, as TADs between different cell types are relative conserved ([Schmitt AD, 2016](#schmitt-ad-2016)), one can use the H1 cell line in case the TAD of corresponding cell type is not available. 

## Epigenetic intensity data
*Irene* requires user to provide [**BigWig**](https://genome.ucsc.edu/goldenpath/help/bigWig.html) format to represent the sequencing density of BS-Seq or ChIP-Seq. User need to create a **data.frame** to indicate the location of the BigWig files, as well as groups and experiment types (as *dataset* column). 
<details><summary>Here is a sample as follows:</summary>

|    | file                                                                                          | group | dataset |
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

</details>

The **BigWig** files for our test cases were converted from wig files retrieved from [NIH Roadmap Epigenomics Project data gateway](https://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/), BLUEPRINT, and [CEEHRC](http://www.epigenomes.ca/), respectively.  Precompiled R objects in our test cases can be retrieved from [here](#precompiled) 

# USAGE
## R environment for the test cases.
The following options were set up during performing our tests. 
```r
options(stringsAsFactors = FALSE)
```

## Load pre-defined regions
The pre-defined promoters and enhancers regions with corresponding IDs (precompiled [hg19](https://github.com/qwang-big/irene-data/blob/master/promenh.hg19.bed), [hg38](https://github.com/qwang-big/irene-data/blob/master/promenh.hg38.bed)) need to be loaded first with: 
```r
bed <- read.table('https://raw.githubusercontent.com/qwang-big/irene-data/master/promenh.hg19.bed')
```

## Import sequencing density data from BigWig files
Read in the table containing file location, group, dataset information as a **data.frame** (named *meta* in the following case), the use 
```r
data <- importBW(meta, bed)
```
to import the density for each regions. *ImportBW* use an external C function from [libBigWig](https://github.com/dpryan79/libBigWig), which is compiled with *Irene* during installation. Its output is equivalent to using *bigWigAverageOverBed* if BigWig files were processed separately, which can be loaded with another [procedure](#bigwigaverageoverbed) instead. 

## Filter out unspecific enhancers
We only took the regions which are more likely to be true enhancers, therefore we use the following function to get the indices which overlapped with enhancer histone marks (H3K4me1 in the following case). The peaks identified by Roadmap Epigenetics Project were retrieved from [http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated](http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated), and run:
```r
i <- filterPeak(c("E003-H3K4me1.narrowPeak","E006-H3K4me1.narrowPeak"), bed, group=c(1,2))
```
 <a id="precompiled">**Precompiled R data files in our test cases (see following) contain all the above mentioned objects (_meta, bed, data, i_) for the following tests**.

| Cancer / primary cells | Controls |
|------------------------|----------|
| Chronic Lymphocytic Leukemia ([CLL](https://github.com/qwang-big/irene-data/blob/master/CLL.hg19.rda)) | Naive B cell |
| Acute Lymphoblastic Leukaemia ([ALL](https://github.com/qwang-big/irene-data/blob/master/ALL.hg38.rda)) | Naive B cell |
| Acute Myeloid Leukaemia ([nkAML](https://github.com/qwang-big/irene-data/blob/master/nkAML.hg38.rda)) | Naive B cell |
| Multiple Myeloma ([MM](https://github.com/qwang-big/irene-data/blob/master/MM.hg38.rda)) | Naive B cell |
| Mantle Cell Lymphoma ([MCL](https://github.com/qwang-big/irene-data/blob/master/MCL.hg38.rda)) | Naive B cell |
| Chronic Lymphocytic Leukemia (mutated) ([mCLL](https://github.com/qwang-big/irene-data/blob/master/mCLL.hg38.rda)) | Naive B cell |
| Colorectal Cancer ([CRC](https://github.com/qwang-big/irene-data/blob/master/CRC.hg19.rda)) | Sigmoid colon |
| Lower Grade Glioma ([LGG](https://github.com/qwang-big/irene-data/blob/master/GLM.hg19.rda)) | Normal brain |
| Papillary Thyroid Cancer ([PTC](https://github.com/qwang-big/irene-data/blob/master/PTC.hg19.rda)) | Normal thyroid |
| Mesenchymal Stem Cells ([MSC](https://github.com/qwang-big/irene-data/blob/master/MSC.hg19.rda)) | Embryonic Stem Cells |
| Neural Progenitor Cells ([NPC](https://github.com/qwang-big/irene-data/blob/master/NPC.hg19.rda)) | Embryonic Stem Cells |
| Trophoblast Stem Cells ([TSC](https://github.com/qwang-big/irene-data/blob/master/TSC.hg19.rda)) | Embryonic Stem Cells |
| H1 BMP4 derived Mesendoderm ([MES](https://github.com/qwang-big/irene-data/blob/master/MES.hg19.rda)) | Embryonic Stem Cells |

*For CLL test case, which is already incorporated in the package, one can simply load necessary dataset with*: 
```r
data(CLL)
```

## Measure combinatorial effect of epigenetic alterations
We use [dPCA](http://www.biostat.jhsph.edu/~hji/dpca/)([Ji H, 2013](#ji-h-2013)) to measure combinatorial effect of epigenetic alterations. The software is already integrated into *Irene* as an external C function. User can select a subset of datasets to study, for the CLL test case:
```r
j <- c(1,2,6,8)
```

Use the following command to read the data which were selected with index *i*: 
```r
res <- dPCA(meta, bed[i,], data[i,], datasets=j, transform=j, normlen=j, verbose=TRUE)
```

or the following if all the genomic regions are to be tested.
```r
res <- dPCA(meta, bed, data, datasets=j, transform=j, normlen=j, verbose=TRUE)
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

The table is sorted descending by PC1

* **Dobs**
The D matrix, which contains the observed differences between the two conditions. This is the data analyzed by dPCA.

* **proj**
Estimated beta coefficients for each dPC.

## Infer promoter-enhancer interaction probabilities according to distances
Use the following function to convert the promoter-enhancer interactions in a given range of TSS (1Mb in the data provided) to probabilities of interaction.

* Load P-E interactions in H1ESC cell lines, the three columns in *H1* data.frame are enhancer IDs, promoter IDs, P-E distances, respectively. 
```r
H1 <- read.table('https://raw.githubusercontent.com/qwang-big/irene-data/master/PEdistances/H1.hg19.pair')
```

* Transform the P-E interactions from bp to Mb:
```r
H1[,3] <- abs(H1[,3]/1000000)
```

* Apply a power-decay function to represent the likelihood of P-E interactions. The choice of power-coefficient is explained in Figure 3.11 and Figure 3.14b of Qi Wang's dissertation. We recommand using -20 for all the test cases. 
```r
H1[,3] <- exp(-20*H1[,3]+1)
```

## Using PageRank to rank the epigenetic alterations from both enhancers and promoters.
The PageRank will take in the alteration contribution from enhancers, resulting the ranks of associated promoters get promoted if targeted by highly altered enhancers. 
```r
res$pg <- pageRank(res$gr, H1)
```

## Get the rank of only promoters as well, which will be used to compare with the PageRank rankings. 
Use the following function to get the order of promoters sorted by PC1, the ones on top of the list indicates higher alteration of the corresponding gene. 
```r
res$pg$prom <- getPromId(res$gr, pc="PC1")
```

## Use literature-derived marker genes as a metric for evaluating the ranking
Cancer and cell-type specific marker genes are listed in Table S2 & 3 of Qi Wang's [disseration](https://github.com/qwang-big/irene-web/blob/master/docs/chapter3.pdf), which are already compiled as an R object and can be loaded with: 
```r
data(markers)
```

The Empirical Cumulative Distribution Function (ECDF) of the marker gene positions in each ranking list can be plotted with:
```r
plotRank(res$pg, markers$CLL)
```
, where we use the CLL marker genes for this test, and the area under the curve (AUC) of each ranking list is described in the legend. 

## Network analysis of the enriched pathways of significant epigenetic alterated genes
Network analysis groups highly-ranked genes according to known gene interaction database, and further searches clusters of genes in biological function databases. Here we loaded Human Protein Reference Database ([HPRD](http://www.hprd.org/)) for grouping genes: 
```r
data(hprd)
```

, and set the edge weights in accordance with the average rank of the two connected genes, which is done by: 
```r
g <- edgeRank(res$pg[[1]], hprd)
```
, where we created a weighted HPRD networks using the edge weights calculated from the PC1 of the genes. 

Afterwards, we create a list of top-ranked sub-networks with desired number, e.g., 15 sub-networks: 
```r
res$gs <- exportMultinets(g, 15)
```

, and the sub-networks were searched in WikiPathways and KEGG for enrichment of biological functions: 
```r
res$ga <- annotNets(res$gs)
```

# Outputs
The following commands need to be executed to generate HTML outputs for interactive exploration of the enriched networks, epigenetic tracks, rank comparisons and benchmarking (ROC curves). 
To begin with, a prefix variable should be set corresponding to the experiment first, as for this test case we set: 
```r
prefix = "CLL"
```

The files are created in the current working directory if not mentioned explicitly. Use the following commands to create a specific folder for the outputs (e.g., a sub-directory in the home folder):
```r
dir.create("~/CLLoutput")
setwd("~/CLLoutput")
```

Write the two ranking list to a CSV file for comparison, significance of the genes from the highest to the lowest are ordered from the top to the bottom of the list.
```r
writeRank(res$pg[[1]], res$pg$prom, prefix)
```

Write the normalized epigenome signals of the two groups, labels of the two groups are provided by user. 
```r
writeData(res$gr, c("CLL", "Bcell"), prefix)
```

, and the sub-networks list was saved to a JSON with each network structure: 
```r
exportJSONnets(res$gs, prefix)
```

In addition, the enriched pathways of the sub-networks are also saved as JSON file:
```r
exportJSONpathways(res$ga, prefix, n=15)
```

Finally, create an index page in the same directory with above JSON files for visualizing in an internet browser:
```r
exportApps(prefix)
```

# Interpreting the results
The above mentioned test cases are hosted in another repository [http://qwang-big.github.io/irene-web](http://qwang-big.github.io/irene-web), in which the outputs are presented in the following four sections: 

## Rank list
For every test case, Irene produced two rank lists: one is from PageRank score (generated from promoter and enhancer scores, as well as P-E interactions, therefore abbreviated as PromEnh), another is from the PC1 scores of promoters (abbreviated as PromOnly). Both lists are contained in a CSV file for downloading, which looks like:

| | PromEnh | PromOnly |
|--|-------------|---------------|
| 1| MED13L | CBWD5 |
| 2| CBLB | CBWD3 |
| 3| VAV3 | CBWD7 |
| 4| LPP | SERF1B |
| 5| NOTCH2 | SERF1A |
| 6| MKLN1 | HIST2H4A |
| 7| ARHGAP15 | HIST2H4B |
| 8| BACH2 | GTF2H2 |
| 9| MGAT5 | NOTCH4 |
...

, where the genes are ordered from the most significant to the least significant in both lists, and the line numbers correspond to the positions in the lists. 

Plotting the positions of each gene in a two-dimensional [graph](http://qwang-big.github.io/irene-web/rank.html?sample=CLL) allows one to inspect the level of significance regarding the promoters or enhancers of a gene, whereas the genes under high enhancer/low promoter regulation are placed at the bottom-right corner of the graph. In the above *CLL* test case, it implies many genes are ranked higher by considering the alteration of enhancers, and the baseline became distorted for those genes who did not receive contribution from enhancers. 
The selected marker genes in benchmarking, as well as from other collections ([COSMIC](https://cancer.sanger.ac.uk/cosmic), [MalaCards](https://www.malacards.org/), [IntOGen](https://www.intogen.org/)) are listed besides the graph, so that the genes are highlighted in red upon clicking the corresponding item. Additional information are shown as mouse-over tooltips or links in the pop-ups. 

## Network enrichment
In these network [presentations](http://qwang-big.github.io/irene-web/net.html?sample=CLL), the enriched sub-networks are ordered based on the average significance of the their nodes, wherein the left-side of the node represents its rank in the *PromEnh* list, and the right-side of the node represents its rank in the *PromOnly* list. Pathways of each sub-network are tested with *EnrichR*, and the ones which are significantly enriched in KEGG and WikiPathways are listed. Upon clicking an item in the drop list of pathways, the corresponding genes in that pathway are highlighted in red. Additional information are shown as mouse-over tooltips or links in the pop-ups. 

## Epigenome heatmap
The epigenome [heatmap](http://qwang-big.github.io/irene-web/browse.html?sample=CLL) serves as an interactive browser of normalized epigenetic signals, wherein tracks from different samples are aligned in accordance with biological conditions and epigenetic marks. In this presentation, the whole genome is divided into bins with fixed length (2kb in this test case). To preserve storage space, only the bins overlaped with the provided promoter/enhancer regions are rendered, leaving the rest of the genome in background color. The top-most track represent the regulatory genomic elements (promoters/enhancers), wherein the colors are corresponding to their PC1 scores and the promoter of the gene name from user input is centered. The triangles above the regulatory elements represent the referred targets from the available resources in the drop list (*GeneHancer* interactions by default). The other tracks represent the data intensities, wherein the colors are corresponding to the normalized epigenetic signals. Additional information are shown as mouse-over tooltips or links in the pop-ups. 

## Hallmark ROC
The receiver operating characteristic ([ROC](http://qwang-big.github.io/irene-web/roc.html?sample=CLL)) curves represent the ECDF of marker genes (curated for benchmarking, [COSMIC](https://cancer.sanger.ac.uk/cosmic), [MalaCards](https://www.malacards.org/), [IntOGen](https://www.intogen.org/)) of the two rank lists (PromEnh and PromOnly). The curves are updated upon clicking the item listed besides the graph. 

# Case studies
Epigenetic and expression data were downloaded from [**EdaccData Release-9**](http://genboree.org/EdaccData/Release-9/experiment-sample/), [**CEEHRC**](http://www.epigenomes.ca/), [**BLUEPRINT**](http://www.blueprint-epigenome.eu/). 

# Acknowledgements
The results presented here are in part based upon data generated by The Canadian Epigenetics, Epigenomics, Environment and Health Research Consortium (CEEHRC) initiative funded by the Canadian Institutes of Health Research (CIHR), Genome BC, and Genome Quebec. Information about CEEHRC and the participating investigators and institutions can be found at [http://www.cihr-irsc.gc.ca/e/43734.html](http://www.cihr-irsc.gc.ca/e/43734.html). 

# References
- <a id="schmitt-ad-2016"></a> Schmitt AD, Hu M, Jung I, et al. A Compendium of Chromatin Contact Maps Reveals Spatially Active Regions in the Human Genome. Cell Reports. 2016;17(8):2042–2059.
- <a id="ji-h-2013"></a> Ji H, Li X, Wang Qf, et al. Differential principal component analysis of ChIP-seq.
Proceedings of the National Academy of Sciences. 2013;110(17):6789–6794.

# License
----
MIT

# Supplymentary Info
- <a id="bigwigaverageoverbed"></a> Constructing **data** object from BigWig files with *bigWigAverageOverBed*: 

# Session Info
```r
pander(sessionInfo(), compact = FALSE)

```

**R version 3.2.2 (2015-08-14)**

**Platform:** x86_64-pc-linux-gnu (64-bit) 

**locale:**
_LC_CTYPE=en_US.UTF-8_, _LC_NUMERIC=C_, _LC_TIME=en_US.UTF-8_, _LC_COLLATE=en_US.UTF-8_, _LC_MONETARY=en_US.UTF-8_, _LC_MESSAGES=en_US.UTF-8_, _LC_PAPER=en_US.UTF-8_, _LC_NAME=C_, _LC_ADDRESS=C_, _LC_TELEPHONE=C_, _LC_MEASUREMENT=en_US.UTF-8_ and _LC_IDENTIFICATION=C_

**attached base packages:** 

* stats4 
* parallel 
* stats 
* graphics 
* grDevices 
* utils 
* datasets 
* methods 
* base 


**other attached packages:** 

* pander(v.0.6.2) 
* irene(v.1.0) 
* enrichR(v.1.0) 
* stringr(v.1.3.1) 
* igraph(v.1.2.2) 
* GenomicRanges(v.1.22.4) 
* GenomeInfoDb(v.1.6.3) 
* IRanges(v.2.4.8) 
* S4Vectors(v.0.8.11) 
* BiocGenerics(v.0.16.1) 


**loaded via a namespace (and not attached):** 

* Rcpp(v.0.12.5) 
* digest(v.0.6.15) 
* R6(v.2.2.2) 
* magrittr(v.1.5) 
* httr(v.1.3.1) 
* stringi(v.1.2.4) 
* zlibbioc(v.1.16.0) 
* XVector(v.0.10.0) 
* preprocessCore(v.1.32.0) 
* rjson(v.0.2.20) 
* tools(v.3.2.2) 
* pkgconfig(v.2.0.2) 
* tcltk(v.3.2.2) 

