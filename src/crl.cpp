#include <Rcpp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
using namespace Rcpp ;
#include <string>
#include <vector>

#include "StringLib.h"
#include "HTSequencingLib.h"

#include "bigWig.h"

//[[Rcpp::export]]
double read_bigwig_impl(std::string bwfname, std::string chrom, int start, int end) {

  //http://stackoverflow.com/questions/347949/how-to-convert-a-stdstring-to-const-char-or-char
  std::vector<char> bwfile(bwfname.begin(), bwfname.end()) ;
  bwfile.push_back('\0') ;

  const char mode = 'r' ;

  bigWigFile_t *bwf = NULL;
  // XXX change NULL to a CURL callback. see libBigWig demos
  bwf = bwOpen(&bwfile[0], NULL, &mode) ;

  if (!bwf)
    exit(1); //stop("Failed to open file: %s\n", bwfname) ;

  double val = 0.0;
  bwOverlappingIntervals_t *intervals = NULL ;

  int nchrom = bwf->cl->nKeys ;
  for (int nc = 0; nc<nchrom; ++nc) {

    char* cur_chrom = bwf->cl->chrom[nc] ;
    std::string cur_chrom_c = cur_chrom ;

    if (!chrom.empty() && chrom != cur_chrom_c) continue ;

    // set maximum boundaries of start / end are not specified
    if (start == -1) start = 0 ;
    if (end == -1) end = bwf->hdr->nBasesCovered ;

    intervals = bwGetValues(bwf, cur_chrom, start, end, 0) ;

    if (!intervals)
      exit(1); //stop("Failed to retreived intervals for %s\n", chrom) ;

    int nint = intervals->l ;
    for(int i=0; i<nint; ++i) {
      val += intervals->value[i] ;
    }

    bwDestroyOverlappingIntervals(intervals) ;
  }

  bwClose(bwf);
  bwCleanup();
  return val ;
}

//[[Rcpp::export]]
List dPCA_main_impl(int nGroupNum, int nDatasetNum, int nSampleNum, int nPaired, int nLociNum, 
		   IntegerVector groupId, IntegerVector datasetId, NumericMatrix repNum, 
		   StringVector sampleName, StringVector lociChr, IntegerVector lociStart,
		   IntegerVector lociEnd, NumericMatrix vdata, 
		  int nTransform, int nColMeanCent, int nColStand, 
		  int nMColMeanCent, int nMColStand, double dSNRCut,
		  int nUsedPCAZ, int nUseRB, double dPeakFDRCut)
{
	/* define */
	struct tagString **vLociChr = NULL; //(struct tagString **)calloc(nLociNum, sizeof(struct tagString *));
	int *vLociStart = lociStart.begin();
	int *vLociEnd = lociEnd.begin();
	
	int *vGroupId = groupId.begin();
	int *vDatasetId = datasetId.begin();
	struct tagString **vSampleName = NULL; //(struct tagString **)calloc(nSampleNum, sizeof(struct tagString *));
	
	int **vRepNum = (int **)calloc(nGroupNum, sizeof(int *));
	float **vData = (float **)calloc(nSampleNum, sizeof(float *));
	double **vM = NULL;
	double dS2;

	/* absolute binding information */
	float **vPeakProb = NULL;
	double **vX1 = NULL;
	double **vX2 = NULL;
	double **vA = NULL;
	
	double **vU = NULL;
	double **vV = NULL;
	double *vH = NULL;
	double *vP = NULL;
	double *vSNR = NULL;
	double **vT = NULL;
	double **vPval = NULL;
	double *vPi = NULL;
	long **vSortId = NULL;
	int nRank;
	double **vRB = NULL;
	
	if(vRepNum == NULL)
	{
		printf("Error: dPCA_Initialize, cannot create replicate number matrix!\n");
		exit(EXIT_FAILURE);
	}
	for(int ni=0; ni< nGroupNum; ni++)
	{
		vRepNum[ni] = (int *)calloc(nDatasetNum, sizeof(int));
		if(vRepNum[ni] == NULL)
		{
			printf("Error: dPCA_Initialize, cannot create replicate number matrix!\n");
			exit(EXIT_FAILURE);
		}
	}
	if(vData == NULL)
	{
		printf("Error: dPCA_Initialize, cannot create data matrix!\n");
		exit(EXIT_FAILURE);
	}
	for(int ni=0; ni<nSampleNum; ni++)
	{
		vData[ni] = (float *)calloc(nLociNum, sizeof(float));
		if(vData[ni] == NULL)
		{
			printf("Error: dPCA_Initialize, cannot create data vector!\n");
			exit(EXIT_FAILURE);
		}
	}

	for(int ni = 0; ni < nGroupNum; ni++) {
	  for(int nj = 0; nj < nDatasetNum; nj++) {
	      vRepNum[ni][nj] = repNum(ni,nj);
	  }
	}
	for(int ni = 0; ni < nSampleNum; ni++) {
	  for(int nj = 0; nj < nLociNum; nj++) {
	      vData[ni][nj] = vdata(nj,ni);
	  }
	}
	//printf("%d\n", vGroupId[nSampleNum-1]);
	//printf("%d\n", vDatasetId[nSampleNum-1]);
	//printf("%d\n", vLociStart[nLociNum-1]);
	//printf("%d\n", vRepNum[nGroupNum-1][nDatasetNum-1]);
	//printf("%f\n", vData[nSampleNum-1][nLociNum-1]);
	
	/* transform, standardization, mean, variance */
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 2: Preprocessing                             */\n");
	printf("/* ------------------------------------------------- */\n");
	if(nPaired == 1)
	{
		dPCA_Preprocess_Paired(vData, nGroupNum, nDatasetNum, nSampleNum, 
			nLociNum, nTransform, nColMeanCent, nColStand, 
			nMColMeanCent, nMColStand, vGroupId, vDatasetId, vRepNum,
			&vM, &dS2, &vX1, &vX2, &vA, nUsedPCAZ, vPeakProb, dPeakFDRCut);
	}
	else
	{
		dPCA_Preprocess_NonPaired(vData, nGroupNum, nDatasetNum, nSampleNum, 
			nLociNum, nTransform, nColMeanCent, nColStand, 
			nMColMeanCent, nMColStand, vGroupId, vDatasetId, vRepNum,
			&vM, &dS2, &vX1, &vX2, &vA, nUsedPCAZ, vPeakProb, dPeakFDRCut);
	}

	/* pca */
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 3: PCA                                       */\n");
	printf("/* ------------------------------------------------- */\n");
	dPCA_SVD(vM, dS2, nGroupNum, nDatasetNum, nLociNum, vRepNum, 
		&vU, &vV, &vH, &vP, &vSNR, dSNRCut, &nRank, &vT, &vPval, nPaired);

	if( (nUsedPCAZ == 1) || (nUseRB == 1) )
	{
		dPCA_RB(nDatasetNum, nLociNum, vM, vV, vPeakProb, dPeakFDRCut, &vRB);
	}

	/* FDR */
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 4: FDR                                       */\n");
	printf("/* ------------------------------------------------- */\n");
	dPCA_FDR(nDatasetNum, nLociNum, nRank, vT, vPval, &vPi, &vSortId);

	NumericVector eigenvalue(nDatasetNum), percent_var(nDatasetNum), SNR(nDatasetNum), sig_prc(nDatasetNum), pc(nLociNum*nDatasetNum), dobs(nLociNum*nDatasetNum);
	for(int ni=0; ni<nDatasetNum; ni++)
		eigenvalue[ni] = vH[ni];
	for(int ni=0; ni<nDatasetNum; ni++)
		percent_var[ni] = vP[ni];
	for(int ni=0; ni<nDatasetNum; ni++)
		SNR[ni] = vSNR[ni];
	for(int ni=0; ni<nDatasetNum; ni++)
		sig_prc[ni] = 1.0-vPi[ni];
	for(int ni=0; ni<nLociNum; ni++)
	{
		for(int nj=0; nj<nDatasetNum; nj++)
		{
			pc[ni*nDatasetNum+nj]   = vU[nj][ni];
			dobs[ni*nDatasetNum+nj] = vM[nj][ni];
		}
	}
	
	/* release memory */
	for(int ni=0; ni<nSampleNum; ni++)
	{
		free(vData[ni]);
		vData[ni] = NULL;

	//	DeleteString(vSampleName[ni]);
	//	vSampleName[ni] = NULL;
	}
	free(vData);
	//free(vSampleName);

	//free(vGroupId);
	//free(vDatasetId);

	for(int ni=0; ni<nGroupNum; ni++)
	{
		free(vRepNum[ni]);
		vRepNum[ni] = NULL;
	}
	free(vRepNum);

	for(int ni=0; ni<nLociNum; ni++)
	{
	//	DeleteString(vLociChr[ni]);
	//	vLociChr[ni] = NULL;
	}
	//free(vLociChr);
	//free(vLociStart);
	//free(vLociEnd);

	for(int ni=0; ni<nDatasetNum; ni++)
	{
		free(vX1[ni]);
		vX1[ni] = NULL;

		free(vX2[ni]);
		vX2[ni] = NULL;

		free(vA[ni]);
		vA[ni] = NULL;

		free(vM[ni]);
		vM[ni] = NULL;

		free(vU[ni]);
		vU[ni] = NULL;

		free(vV[ni]);
		vV[ni] = NULL;

		free(vT[ni]);
		vT[ni] = NULL;

		free(vPval[ni]);
		vPval[ni] = NULL;

		free(vSortId[ni]);
		vSortId[ni] = NULL;
	}
	free(vX1);
	free(vX2);
	free(vA);
	free(vM);
	free(vU);
	free(vV);
	free(vH);
	free(vP);
	free(vSNR);
	free(vT);
	free(vPval);
	free(vPi);
	free(vSortId);

	if( (nUsedPCAZ == 1) || (nUseRB == 1) )
	{
		for(int ni=0; ni<nDatasetNum; ni++)
		{
			free(vPeakProb[ni]);
			vPeakProb[ni] = NULL;

			free(vRB[ni]);
			vRB[ni] = NULL;
		}
		free(vPeakProb);
		free(vRB);
	}

	/* return */
	return Rcpp::List::create(Rcpp::Named("eigenvalue") = eigenvalue,Rcpp::Named("percent_var") = percent_var,Rcpp::Named("SNR") = SNR,
                          Rcpp::Named("sig_prc") = sig_prc,Rcpp::Named("s2") = dS2,Rcpp::Named("PC") = pc,Rcpp::Named("Dobs") = dobs);
}
