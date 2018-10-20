/* ----------------------------------------------------------------------- */
/*  HTSequencingLib.c : implementation of the high throughput sequencing   */
/*  library                                                                */
/*  Author : Ji HongKai ; Time: 2007.10                                    */
/* ----------------------------------------------------------------------- */

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "limits.h"

#include "MathLib.h"
#include "MatrixLib.h"
#include "RandomLib.h"
#include "StringLib.h"
#include "GenomeLib.h"
#include "AffyLib.h"
#include "TilingArrayLib.h"
#include "HTSequencingLib.h"
#include "CisChromeLib.h"


/* ----------------------------------------------------------------------- */ 
/*  CisChrome_Main()                                                       */
/*  CisChrome: analyzing chromatin signals at DNA motif sites.             */
/* ----------------------------------------------------------------------- */ 
int CisChrome_Main(char strParamFile[])
{
	/* define */
	char strWorkFolder[MED_LINE_LENGTH];
	char strProjectTitle[MED_LINE_LENGTH];
	char strMotifFolder[MED_LINE_LENGTH];
	char strMotifLibrary[MED_LINE_LENGTH];
	int nMotifNum;
	struct tagString **vTFName;
	struct tagRefGene **vTFRefGene;
	struct tagString **vMotifFilePath;
	int *vTFEntrezID;
	char strSpecies[LINE_LENGTH];

	int nGridNum = 10;

	int nGroupNum = 0;
	int nDatasetNum = 0;
	int nSampleNum = 0;
	int nPaired = 0;
	int nBinSize = 100;
	int nWinSize = 500;
	int nTwosided = 1;
	int *vExtLen = NULL;
	int *vGroupId = NULL;
	int *vDatasetId = NULL;
	struct tagString **vSampleName;
	struct tagString **vSampleFilePath;
	int **vRepNum = NULL;
	int ni;

	int nChrNum = 0;
	struct tagString **vChrName = NULL;
	struct INTMATRIX *pChrLen = NULL;

	float **vMotifStat = NULL;
	float **vMotifCDF = NULL;
	double *vMinT,*vMaxT;
	double dFDRCut = 0.05;
	int *vTotalSiteNum;
	int *vActiveSiteNum;
	double *vActiveSitePrc;
	double *vTFExpress;
	int nUp = 2000;
	int nDown = 1000;
	int nTSS = 1;

	/* ---------------------------------------------------------*/
	/* STEP1: Read Parameters.                                  */
	/* ---------------------------------------------------------*/
	CisChrome_ReadParameters(strParamFile, strWorkFolder, strProjectTitle, 
		strMotifFolder, strMotifLibrary, &nGroupNum, &nDatasetNum, &nSampleNum, 
		&nPaired, &nBinSize, &nWinSize,	&nTwosided, &vSampleName, 
		&vSampleFilePath, &vGroupId, &vDatasetId, &vExtLen, &vRepNum,
		&nGridNum, &dFDRCut, &nTSS, &nUp, &nDown);

	CisChrome_ReadTFInfo(strMotifLibrary, &nMotifNum, &vTFName, 
		&vTFRefGene, &vMotifFilePath, &vTFEntrezID, strSpecies);

	/* ---------------------------------------------------------*/
	/* STEP2: Initialize.                                       */
	/* ---------------------------------------------------------*/
	vMinT = NULL;
	vMinT = (double *)calloc(nDatasetNum, sizeof(double));
	if(vMinT == NULL)
	{
		printf("Error: CisChrome_Main, cannot allocate memory for motif minT statistics!\n");
		exit(EXIT_FAILURE);
	}
	vMaxT = NULL;
	vMaxT = (double *)calloc(nDatasetNum, sizeof(double));
	if(vMaxT == NULL)
	{
		printf("Error: CisChrome_Main, cannot allocate memory for motif maxT statistics!\n");
		exit(EXIT_FAILURE);
	}
	CisChrome_GenomeProfile(nGroupNum, nDatasetNum, nSampleNum, nPaired, nTwosided,
		vGroupId, vDatasetId, vExtLen, vRepNum,
		vSampleName, vSampleFilePath,
		nBinSize, nWinSize, strWorkFolder, strProjectTitle,
		&nChrNum, &vChrName, &pChrLen, nGridNum, vMinT, vMaxT);

	/* ---------------------------------------------------------*/
	/* STEP3: Process motifs one-by-one.                        */
	/* ---------------------------------------------------------*/
	vMotifStat = NULL;
	vMotifStat = (float **)calloc(nMotifNum, sizeof(float *));
	if(vMotifStat == NULL)
	{
		printf("Error: CisChrome_Main, cannot allocate memory for motif statistics!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nMotifNum; ni++)
	{
		vMotifStat[ni] = (float *)calloc(nGridNum*nDatasetNum, sizeof(float));
		if(vMotifStat[ni] == NULL)
		{
			printf("Error: CisChrome_Main, cannot allocate memory for motif statistics!\n");
			exit(EXIT_FAILURE);
		}
	}
	vMotifCDF = NULL;
	vMotifCDF = (float **)calloc(nMotifNum, sizeof(float *));
	if(vMotifCDF == NULL)
	{
		printf("Error: CisChrome_Main, cannot allocate memory for motif CDF statistics!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nMotifNum; ni++)
	{
		vMotifCDF[ni] = (float *)calloc(nGridNum*nDatasetNum, sizeof(float));
		if(vMotifCDF[ni] == NULL)
		{
			printf("Error: CisChrome_Main, cannot allocate memory for motif statistics!\n");
			exit(EXIT_FAILURE);
		}
	}

	vTotalSiteNum = NULL;
	vTotalSiteNum = (int *)calloc(nMotifNum, sizeof(int));
	if(vTotalSiteNum == NULL)
	{
		printf("Error: CisChrome_Main, cannot allocate memory for total motif site number!\n");
		exit(EXIT_FAILURE);
	}

	vActiveSiteNum = NULL;
	vActiveSiteNum = (int *)calloc(nMotifNum*nDatasetNum, sizeof(int));
	if(vActiveSiteNum == NULL)
	{
		printf("Error: CisChrome_Main, cannot allocate memory for active motif site number!\n");
		exit(EXIT_FAILURE);
	}

	vActiveSitePrc = NULL;
	vActiveSitePrc = (double *)calloc(nMotifNum*nDatasetNum, sizeof(double));
	if(vActiveSitePrc == NULL)
	{
		printf("Error: CisChrome_Main, cannot allocate memory for active motif site percentage!\n");
		exit(EXIT_FAILURE);
	}

	vTFExpress = NULL;
	vTFExpress = (double *)calloc(nMotifNum*nDatasetNum, sizeof(double));
	if(vTFExpress == NULL)
	{
		printf("Error: CisChrome_Main, cannot allocate memory for TF expression statistics!\n");
		exit(EXIT_FAILURE);
	}

	CisChrome_FindActiveTFBS(strWorkFolder, strProjectTitle,
		nDatasetNum, nBinSize, nWinSize, nTwosided,
		nChrNum, vChrName, pChrLen, nGridNum,
		strMotifFolder, nMotifNum, vTFName, vMotifFilePath, 
		vMotifStat, vMotifCDF, vMinT, vMaxT, dFDRCut,
		vTotalSiteNum, vActiveSiteNum, vActiveSitePrc,
		vTFExpress, vTFRefGene, nUp, nDown, nTSS);

	/* ---------------------------------------------------------*/
	/* STEP4: Export results.                                   */
	/* ---------------------------------------------------------*/
	CisChrome_Export(strWorkFolder, strProjectTitle, nMotifNum, nDatasetNum,
		nGridNum, vTFName, vTotalSiteNum, vActiveSiteNum, vActiveSitePrc,
		vMotifStat, vMinT, vMaxT, vTFExpress, 0);

	CisChrome_Export(strWorkFolder, strProjectTitle, nMotifNum, nDatasetNum,
		nGridNum, vTFName, vTotalSiteNum, vActiveSiteNum, vActiveSitePrc,
		vMotifCDF, vMinT, vMaxT, vTFExpress, 1);

	/* ---------------------------------------------------------*/
	/* STEP5: Release memory.                                   */
	/* ---------------------------------------------------------*/
	for(ni=0; ni<nSampleNum; ni++)
	{
		DeleteString(vSampleName[ni]);
		vSampleName[ni] = NULL;

		DeleteString(vSampleFilePath[ni]);
		vSampleFilePath[ni] = NULL;
	}
	free(vSampleName);
	free(vSampleFilePath);

	free(vExtLen);
	free(vGroupId);
	free(vDatasetId);

	for(ni=0; ni<nGroupNum; ni++)
	{
		free(vRepNum[ni]);
		vRepNum[ni] = NULL;
	}
	free(vRepNum);

	for(ni=0; ni<nChrNum; ni++)
	{
		DeleteString(vChrName[ni]);
		vChrName[ni] = NULL;
	}
	free(vChrName);
	DestroyIntMatrix(pChrLen);

	for(ni=0; ni<nMotifNum; ni++)
	{
		DeleteString(vTFName[ni]);
		vTFName[ni] = NULL;

		RefGeneDestroy(vTFRefGene[ni]);
		vTFRefGene[ni] = NULL;

		DeleteString(vMotifFilePath[ni]);
		vMotifFilePath[ni] = NULL;

		free(vMotifStat[ni]);
		vMotifStat[ni] = NULL;

		free(vMotifCDF[ni]);
		vMotifCDF[ni] = NULL;
	}
	free(vTFName);
	free(vTFRefGene);
	free(vMotifFilePath);
	free(vTFEntrezID);
	free(vMotifStat);
	free(vMotifCDF);
	free(vMinT);
	free(vMaxT);

	free(vTotalSiteNum);
	free(vActiveSiteNum);
	free(vActiveSitePrc);
	free(vTFExpress);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  CisChrome_ReadParameters()                                             */
/*  Read parameters.                                                       */
/* ----------------------------------------------------------------------- */
int CisChrome_ReadParameters(char strParamFile[], char *strWorkFolder, 
		char *strProjectTitle, char *strMotifFolder, char *strMotifLibrary, 
		int *pGroupNum, int *pDatasetNum, int *pSampleNum, int *pPaired, 
		int *pBinSize, int *pWinSize, int *pTwosided, 
		struct tagString ***pvSampleName, struct tagString ***pvSampleFilePath, 
		int **pvGroupId, int **pvDatasetId, int **pvExtLen, int ***pvRepNum,
		int *pGridNum, double *pFDRCut, int *pTSS, int *pUp, int *pDown)
{
	/* define */
	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	char *chp,*chp2;
	int ni;
	char strSampleName[MED_LINE_LENGTH];
	char strSamplePath[MED_LINE_LENGTH];
	int nGi,nDi;

	/* init */
	*pDatasetNum = 0;
	*pSampleNum = 0;

	/* open */
	fpIn = NULL;
	fpIn = fopen(strParamFile, "rt");
	if(fpIn == NULL)
	{
		printf("Error: CisChrome_ReadParameters, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	/* read */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		if(strstr(strLine, "[WorkingFolder]") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			strcpy(strWorkFolder, chp);
			AdjustDirectoryPath(strWorkFolder);
		}
		else if(strstr(strLine, "[ProjectTitle]") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			strcpy(strProjectTitle, chp);
		}
		else if(strstr(strLine, "[MotifFolder]") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			strcpy(strMotifFolder, chp);
			AdjustDirectoryPath(strMotifFolder);
		}
		else if(strstr(strLine, "[MotifLibrary]") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			strcpy(strMotifLibrary, chp);
		}
		else if(strstr(strLine, "[TFETSS]") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			*pTSS = atoi(chp);
		}
		else if(strstr(strLine, "[TFEUP]") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			*pUp = atoi(chp);
		}
		else if(strstr(strLine, "[TFEDOWN]") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			*pDown = atoi(chp);
		}
		else if(strstr(strLine, "[GridNum]") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			*pGridNum = atoi(chp);
		}
		else if(strstr(strLine, "[GroupNum]") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			*pGroupNum = atoi(chp);

			if(*pGroupNum > 2)
			{
				printf("Error: CisChrome_ReadParameters, currently we do not support more than 2 groups/conditions, sorry!\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strstr(strLine, "[DatasetNum]") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			*pDatasetNum = atoi(chp);
		}
		else if(strstr(strLine, "[SampleNum]") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			*pSampleNum = atoi(chp);
		}
		else if(strstr(strLine, "[Paired]") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			*pPaired = atoi(chp);
		}
		else if(strstr(strLine, "[BinSize]") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			*pBinSize = atoi(chp);
		}
		else if(strstr(strLine, "[WinSize]") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			*pWinSize = atoi(chp);
		}
		else if(strstr(strLine, "[TwoSided]") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			*pTwosided = atoi(chp);
		}
		else if(strstr(strLine, "[FDRCut]") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			*pFDRCut = atof(chp);
		}
		else if(strstr(strLine, "[Samples]") == strLine)
		{
			if(*pSampleNum == 0)
			{
				printf("Warning: CisChrome_ReadParameters, no samples available!\n");
				return PROC_SUCCESS;
			}

			*pvSampleName = NULL;
			*pvSampleName = (struct tagString **)calloc(*pSampleNum, sizeof(struct tagString *));
			if(*pvSampleName == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot create sample name vector!\n");
				exit(EXIT_FAILURE);
			}

			*pvSampleFilePath = NULL;
			*pvSampleFilePath = (struct tagString **)calloc(*pSampleNum, sizeof(struct tagString *));
			if(*pvSampleFilePath == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot create sample file path vector!\n");
				exit(EXIT_FAILURE);
			}

			*pvGroupId = NULL;
			*pvGroupId = (int *)calloc(*pSampleNum, sizeof(int));
			if(*pvGroupId == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot create group id vector!\n");
				exit(EXIT_FAILURE);
			}

			*pvDatasetId = NULL;
			*pvDatasetId = (int *)calloc(*pSampleNum, sizeof(int));
			if(*pvDatasetId == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot create dataset id vector!\n");
				exit(EXIT_FAILURE);
			}

			*pvExtLen = NULL;
			*pvExtLen = (int *)calloc(*pSampleNum, sizeof(int));
			if(*pvExtLen == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot create extension length vector!\n");
				exit(EXIT_FAILURE);
			}

			*pvRepNum = NULL;
			*pvRepNum = (int **)calloc(*pGroupNum, sizeof(int *));
			if(*pvRepNum == NULL)
			{
				printf("Error: CisChrome_ReadParameters, cannot create replicate number matrix!\n");
				exit(EXIT_FAILURE);
			}
			for(ni=0; ni< (*pGroupNum); ni++)
			{
				(*pvRepNum)[ni] = (int *)calloc(*pDatasetNum, sizeof(int));
				if((*pvRepNum)[ni] == NULL)
				{
					printf("Error: CisChrome_ReadParameters, cannot create replicate number matrix!\n");
					exit(EXIT_FAILURE);
				}
			}

			ni = 0;
			while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;
				if(strLine[0] == '#')
					continue;

				chp2 = strchr(strLine, '\t');
				if(chp2 == NULL)
				{
					printf("Error: CisChrome_ReadParameters, parameter file has wrong format, sample information should be tab-delimited!\n");
					exit(EXIT_FAILURE);
				}
				*chp2 = '\0';
				strcpy(strSampleName, strLine);
				StringAddTail((*pvSampleName)+ni, strSampleName);
				chp = chp2+1;

				chp2 = strchr(chp, '\t');
				if(chp2 == NULL)
				{
					printf("Error: CisChrome_ReadParameters, parameter file has wrong format, sample information should be tab-delimited!\n");
					exit(EXIT_FAILURE);
				}
				*chp2 = '\0';
				nGi = atoi(chp);
				(*pvGroupId)[ni] = nGi;
				chp = chp2+1;

				chp2 = strchr(chp, '\t');
				if(chp2 == NULL)
				{
					printf("Error: CisChrome_ReadParameters, parameter file has wrong format, sample information should be tab-delimited!\n");
					exit(EXIT_FAILURE);
				}
				*chp2 = '\0';
				nDi = atoi(chp);
				(*pvDatasetId)[ni] = nDi;
				chp = chp2+1;

				chp2 = strchr(chp, '\t');
				if(chp2 == NULL)
				{
					printf("Error: CisChrome_ReadParameters, parameter file has wrong format, sample information should be tab-delimited!\n");
					exit(EXIT_FAILURE);
				}
				*chp2 = '\0';
				(*pvExtLen)[ni] = atoi(chp);
				chp = chp2+1;

				StrTrimLeft(chp);
				strcpy(strSamplePath, chp);
				StringAddTail((*pvSampleFilePath)+ni, strSamplePath);

				(*pvRepNum)[nGi-1][nDi-1] += 1;

				ni++;
			}

			if(ni != (*pSampleNum))
			{
				printf("Error: CisChrome_ReadParameters, sample number not match!\n");
				exit(EXIT_FAILURE);
			}

		}
	}

	/* close file */
	fclose(fpIn);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_ReadTFInfo()                                                 */
/*  Read TF information.                                                   */
/* ----------------------------------------------------------------------- */
int CisChrome_ReadTFInfo(char strMotifLibrary[], int *pMotifNum, 
	struct tagString ***pvTFName, struct tagRefGene ***pvTFRefGene,
	struct tagString ***pvMotifFilePath, int **pvTFEntrezID, char strSpecies[])
{
	/* define */
	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	char *chp;
	int nMotifId;

	/* init */
	*pMotifNum = 0;
	*pvTFName = NULL;
	*pvTFRefGene = NULL;
	*pvMotifFilePath = NULL;
	*pvTFEntrezID = NULL;

	/* open */
	fpIn = NULL;
	fpIn = fopen(strMotifLibrary, "rt");
	if(fpIn == NULL)
	{
		printf("Error: CisChrome_ReadTFInfo, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	/* read */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		if(strstr(strLine, "SPECIES") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadTFInfo, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			strcpy(strSpecies, chp);
		}
		else if(strstr(strLine, "TF_NUM") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadTFInfo, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			*pMotifNum = atoi(chp);

			nMotifId = -1;

			*pvTFName = NULL;
			*pvTFName = (struct tagString **)calloc(*pMotifNum, sizeof(struct tagString *));
			if( *pvTFName == NULL )
			{
				printf("Error: CisChrome_ReadTFInfo, cannot allocate memory for TF names!\n");
				exit(EXIT_FAILURE);
			}

			*pvTFRefGene = NULL;
			*pvTFRefGene = (struct tagRefGene **)calloc(*pMotifNum, sizeof(struct tagRefGene *));
			if( *pvTFRefGene == NULL )
			{
				printf("Error: CisChrome_ReadTFInfo, cannot allocate memory for TF names!\n");
				exit(EXIT_FAILURE);
			}

			*pvMotifFilePath = NULL;
			*pvMotifFilePath = (struct tagString **)calloc(*pMotifNum, sizeof(struct tagString *));
			if( *pvMotifFilePath == NULL )
			{
				printf("Error: CisChrome_ReadTFInfo, cannot allocate memory for motif site paths!\n");
				exit(EXIT_FAILURE);
			}

			*pvTFEntrezID = NULL;
			*pvTFEntrezID = (int *)calloc(*pMotifNum, sizeof(int));
			if( *pvTFEntrezID == NULL )
			{
				printf("Error: CisChrome_ReadTFInfo, cannot allocate memory for TF Entrez ID!\n");
				exit(EXIT_FAILURE);
			}

			
		}
		else if(strstr(strLine, ">") == strLine)
		{
			nMotifId++;
		}
		else if(strstr(strLine, "NAME") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadTFInfo, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			StringAddTail((*pvTFName)+nMotifId, chp);
		}
		else if(strstr(strLine, "ENTREZ") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadTFInfo, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			(*pvTFEntrezID)[nMotifId] = atoi(chp);
		}
		else if(strstr(strLine, "COORDINATE") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadTFInfo, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			(*pvTFRefGene)[nMotifId] = RefGeneCreate();
			RefFlatInit((*pvTFRefGene)[nMotifId], chp, strSpecies);
		}
		else if(strstr(strLine, "MOTIFSITES") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp == NULL)
			{
				printf("Error: CisChrome_ReadTFInfo, cannot find = separator in the parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			StrTrimLeft(chp);
			StringAddTail((*pvMotifFilePath)+nMotifId, chp);
		}
	}

	/* close file */
	fclose(fpIn);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_GenomeProfile()                                              */
/*  Create Genome-wide enrichment profiles.                                */
/* ----------------------------------------------------------------------- */
int CisChrome_GenomeProfile(int nGroupNum, int nDatasetNum, int nSampleNum, 
		int nPaired, int nTwosided, int *vGroupId, int *vDatasetId, int *vExtLen, 
		int **vRepNum, struct tagString **vSampleName, struct tagString **vSampleFilePath,
		int nBinSize,  int nWinSize, char strWorkFolder[], char strProjectTitle[],
		int *pChrNum, struct tagString ***vChrName, struct INTMATRIX **pChrLen,
		int nGridNum, double *pMinT, double *pMaxT)
{
	/* define */
	struct DOUBLEMATRIX *pFileReadCount = NULL;
	struct DOUBLEMATRIX *pFileNormFactor = NULL;
	int nBaseLineFileId = 0;
	int nBaseReadNum = 0;
	int nExportBAR = 0;
	int nResult;

	/* ---------------------------------------------------------*/
	/* STEP1: First scan, find chromosome number and length     */
	/*        Find read number for each sample                  */
	/* ---------------------------------------------------------*/
	printf("/* ------------------------------------------------- */\n");
	printf("/* Initialize                                        */\n");
	printf("/* ------------------------------------------------- */\n");
	*pChrNum = 0;
	*vChrName = NULL;
	*pChrLen = NULL;
	nResult = CisChrome_GenomeProfile_Initialize(nSampleNum, 
			vSampleName, vSampleFilePath, vExtLen, nBinSize,
			&pFileReadCount, &pFileNormFactor, 
			&nBaseReadNum, &nBaseLineFileId, 
			pChrNum, vChrName, pChrLen);
	
	/* ---------------------------------------------------------*/
	/* STEP2: Count reads for genomic bins                      */
	/* ---------------------------------------------------------*/
	printf("\n");
	printf("/* ------------------------------------------------- */\n");
	printf("/* Count reads for genomic bins                      */\n");
	printf("/* ------------------------------------------------- */\n");
	nResult = CisChrome_GenomeProfile_CountBinReads(nSampleNum, 
		vSampleName, vSampleFilePath, vExtLen, nBinSize,
		pFileNormFactor, strWorkFolder, strProjectTitle,  
		(*pChrNum), (*vChrName), (*pChrLen), nExportBAR);
	

	/* ---------------------------------------------------------*/
	/* STEP3: Compute summary statistic for genomic bins        */
	/* ---------------------------------------------------------*/
	printf("\n");
	printf("/* ------------------------------------------------- */\n");
	printf("/* Compute summary statistics for bins               */\n");
	printf("/* ------------------------------------------------- */\n");
	nResult = CisChrome_ComputeBinStats(nGroupNum, nDatasetNum, nSampleNum, 
		nPaired, nTwosided, vGroupId, vDatasetId, vRepNum, vSampleName, nBinSize, nWinSize,
		strWorkFolder, strProjectTitle, (*pChrNum), (*vChrName), (*pChrLen),
		nGridNum, pMinT, pMaxT);

	/* release memory */
	DestroyDoubleMatrix(pFileReadCount);
	DestroyDoubleMatrix(pFileNormFactor);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_GenomeProfile_Initialize()                                   */
/*  Initialize genomewide enrichment profiles.                             */
/* ----------------------------------------------------------------------- */
int CisChrome_GenomeProfile_Initialize(int nSampleNum, 
		struct tagString **vSampleName, struct tagString **vSampleFilePath, 
		int *vExtLen, int nBinSize, struct DOUBLEMATRIX **ppFileReadCount, 
		struct DOUBLEMATRIX **ppFileNormFactor, 
		int *pBaseReadNum, int *pBaseLineFileId, 
		int *pChrNum, struct tagString ***vvChrName, struct INTMATRIX **ppChrLen)
{
	/* define */
	int nMaxChrNum = 65535;	
	int ni,nId;
	int nChrNum = 0;
	int nReadCount;
	double dBaseCount;
	struct tagProbeGenomeInfo **vChrList = NULL;

	FILE *fpRead;
	char strLine[LONG_LINE_LENGTH];
	struct DOUBLEMATRIX *pSampleReadCount;
	struct DOUBLEMATRIX *pSampleReadCountSort;
	struct LONGMATRIX *pSampleReadCountIdx;
	char strChr[LINE_LENGTH];
	int nPos;
	int nStrand;

	/* Count file number */
	printf("No. of samples = %d\n", nSampleNum);
	
	*ppFileReadCount = NULL;
	*ppFileReadCount = CreateDoubleMatrix(1, nSampleNum);
	if(*ppFileReadCount == NULL)
	{
		printf("Error: CisChrome_GenomeProfile_Initialize, cannot create vector for file read counts!\n");
		exit(EXIT_FAILURE);
	}

	pSampleReadCount = NULL;
	pSampleReadCount = CreateDoubleMatrix(1, nSampleNum);
	if(pSampleReadCount == NULL)
	{
		printf("Error: CisChrome_GenomeProfile_Initialize, cannot create vector for file read counts!\n");
		exit(EXIT_FAILURE);
	}

	*ppFileNormFactor = NULL;
	*ppFileNormFactor = CreateDoubleMatrix(1, nSampleNum);
	if(*ppFileNormFactor == NULL)
	{
		printf("Error: CisChrome_GenomeProfile_Initialize, cannot create vector for file normalizing factor!\n");
		exit(EXIT_FAILURE);
	}

	/* process files one by one */
	vChrList = NULL;
	vChrList = (struct tagProbeGenomeInfo **)calloc(nMaxChrNum, sizeof(struct tagProbeGenomeInfo *));
	if(vChrList == NULL)
	{
		printf("Error: CisChrome_GenomeProfile_Initialize, cannot create vector for chromosomes!\n");
		exit(EXIT_FAILURE);
	}
	nChrNum = 0;

	for(ni=0; ni<nSampleNum; ni++)
	{	
		printf("  Processing %s ...\n", vSampleFilePath[ni]->m_pString);
		nReadCount = 0;

		/* IF ALN FILE */
		fpRead = NULL;
		fpRead = fopen(vSampleFilePath[ni]->m_pString, "r");
		if(fpRead == NULL)
		{
			printf("Error: CisChrome_GenomeProfile_Initialize, cannot open read alignment file!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine, LONG_LINE_LENGTH, fpRead) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			if(strLine[0] == '#')
				continue;

			SeqClust_ReadPos_From_Aln(strLine, strChr, &nPos, &nStrand);

			/* find chromosome and update chromosome length */
			nId = SeqClust_FindChrInfo(vChrList, nMaxChrNum, &nChrNum, strChr);
			if( (nId < 0) || (nId >= nMaxChrNum) )
			{
				printf("Error: CisChrome_GenomeProfile_Initialize, cannot find the matching chromosome!\n");
				exit(EXIT_FAILURE);
			}
			if(nPos > vChrList[nId]->nPos)
			{
				vChrList[nId]->nPos = nPos;
			}
			
			nReadCount++;
		}	

		/* close alignment file */
		fclose(fpRead);

		/* TODO: IF BAM FILE */

		/* update file information */
		pSampleReadCount->pMatElement[ni] = nReadCount;
		(*ppFileReadCount)->pMatElement[ni] = nReadCount;
	}
	
	if( ni != nSampleNum )
	{
		printf("Error: CisChrome_GenomeProfile_Initialize, inconsistent sample number!\n");
		exit(EXIT_FAILURE);
	}

	/* return chromosome name and length */
	*pChrNum = nChrNum;
	*vvChrName = (struct tagString **)calloc(nChrNum, sizeof(struct tagString *));
	if(*vvChrName == NULL)
	{
		printf("Error: CisChrome_GenomeProfile_Initialize, cannot create memory for chromosome name!\n");
		exit(EXIT_FAILURE);
	}

	*ppChrLen = NULL;
	*ppChrLen = CreateIntMatrix(1, nChrNum);
	if(*ppChrLen == NULL)
	{
		printf("Error: CisChrome_GenomeProfile_Initialize, cannot create memory for chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nChrNum; ni++)
	{
		(*vvChrName)[ni] = vChrList[ni]->pProbe;
		vChrList[ni]->pProbe = NULL;
		(*ppChrLen)->pMatElement[ni] = vChrList[ni]->nPos + 1;

		ProbeGenomeInfoDestroy(vChrList+ni);
	}
	free(vChrList);

	if(ni != nChrNum)
	{
		printf("Error: CisChrome_GenomeProfile_Initialize, inconsistent chromosome number!\n");
		exit(EXIT_FAILURE);
	}

	/* compute normalizing factor */
	pSampleReadCountSort = NULL;
	pSampleReadCountIdx = NULL;
	DMSORTMERGEA_0(pSampleReadCount, &pSampleReadCountSort, &pSampleReadCountIdx);
	
	dBaseCount = pSampleReadCountSort->pMatElement[0];
	*pBaseLineFileId = pSampleReadCountIdx->pMatElement[0];
	*pBaseReadNum = (int)dBaseCount;

	for(ni=0; ni<nSampleNum; ni++)
	{
		if( (*ppFileReadCount)->pMatElement[ni] <= 0)
		{
			printf("Error: CisChrome_GenomeProfile_Initialize, sample read count for the %d-th file is zero, cannot perform normalization!\n", ni+1);
			exit(EXIT_FAILURE);
		}
		(*ppFileNormFactor)->pMatElement[ni] = dBaseCount*nBinSize/(*ppFileReadCount)->pMatElement[ni]/vExtLen[ni];
	}


	/* print information */
	printf("\nFile\tNo_of_Reads\tScaling_Factor\tExtension_Length\n");
	for(ni=0; ni<nSampleNum; ni++)
	{
		printf("%s\t%d\t%f\t%d\n", (vSampleName)[ni]->m_pString, (int)((*ppFileReadCount)->pMatElement[ni]), 
			(*ppFileNormFactor)->pMatElement[ni], (vExtLen)[ni]);
	}

	printf("\nNo. of chromosomes = %d\n", *pChrNum);
	printf("\tChromosome\tMax_Coordinate\n");
	for(ni=0; ni<(*pChrNum); ni++)
	{
		printf("\t%s\t%d\n", (*vvChrName)[ni]->m_pString, (*ppChrLen)->pMatElement[ni]);
	}

	/* release memory */
	DestroyDoubleMatrix(pSampleReadCount);
	DestroyDoubleMatrix(pSampleReadCountSort);
	DestroyLongMatrix(pSampleReadCountIdx);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_GenomeProfile_CountBinReads()                                */
/*  Count reads for genomic bins                                           */
/* ----------------------------------------------------------------------- */ 
int CisChrome_GenomeProfile_CountBinReads(int nSampleNum, 
		struct tagString **vSampleName, struct tagString **vSampleFilePath, 
		int *vExtLen, int nBinSize,	struct DOUBLEMATRIX *pFileNormFactor,
		char strWorkFolder[], char strProjectTitle[],  
		int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen,
		int nExportBAR)
{
	/* define */
	int ni;
	int nResult;

	for(ni=0; ni<nSampleNum; ni++)
	{
		printf("  Processing %s ...\n", vSampleFilePath[ni]->m_pString);
		nResult = CisChrome_GenomeProfile_CountBinReads_SingleFile(vSampleFilePath[ni]->m_pString,
			strWorkFolder, vSampleName[ni]->m_pString, nChrNum, vChrName, pChrLen, nBinSize,
			vExtLen[ni], pFileNormFactor->pMatElement[ni], nExportBAR);
	}
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_GenomeProfile_CountBinReads_SingleFile()                     */
/*  Count reads for genomic bins for a single sample file                  */
/* ----------------------------------------------------------------------- */ 
int CisChrome_GenomeProfile_CountBinReads_SingleFile(char strInputFile[], char strOutputPath[],
				char strOutputFile[], int nChrNum, struct tagString **vChrName, 
				struct INTMATRIX *pChrLen, int nBinSize, 		
				int nExtLen, double dNormFactor, int nExportBAR)
{
	/* define */
	float **vBinC = NULL;
	int ni,nk,nx;
	int *vBinNum;
	int nBinNum;
	char strFileName[LINE_LENGTH];
	char strOutFileName[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	FILE *fpRead;
	char strChr[LINE_LENGTH];
	int nPos,nEndPos;
	int nStrand;
	int nChrId,nBinId,nEndBinId;
	float fTemp;

	/* init */
	GetFileName(strInputFile, strFileName);
	vBinC = NULL;
	vBinC = (float **)calloc(nChrNum, sizeof(float *));
	if(vBinC == NULL)
	{
		printf("Error: CisChrome_GenomeProfile_CountBinReads_SingleFile, cannot create vector for genomic bins!\n");
		exit(EXIT_FAILURE);
	}

	vBinNum = NULL;
	vBinNum = (int *)calloc(nChrNum, sizeof(int));
	if(vBinNum == NULL)
	{
		printf("Error: CisChrome_GenomeProfile_CountBinReads_SingleFile, cannot create vector for genomic bin size!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nChrNum; ni++)
	{
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		vBinC[ni] = (float *)calloc(nBinNum, sizeof(float));
		if(vBinC[ni] == NULL)
		{
			printf("Error: CisChrome_GenomeProfile_CountBinReads_SingleFile, insufficient memory for creating genomic bins, try a larger bin size!\n");
			exit(EXIT_FAILURE);
		}		

		vBinNum[ni] = nBinNum;
	}

	/* IF ALN FILE */
	fpRead = NULL;
	fpRead = fopen(strInputFile, "r");
	if(fpRead == NULL)
	{
		printf("Error: CisChrome_GenomeProfile_CountBinReads_SingleFile, cannot open read alignment file!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, LONG_LINE_LENGTH, fpRead) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		SeqClust_ReadPos_From_Aln(strLine, strChr, &nPos, &nStrand);

		/* find chromosome and update chromosome length */
		nChrId = SeqClust_FindChr(vChrName, nChrNum, strChr);
		if( (nChrId < 0) || (nChrId >= nChrNum) )
		{
			printf("Warning: CisChrome_GenomeProfile_CountBinReads_SingleFile, cannot find the matching chromosome for %s:%d!\n",strChr, nPos);
			continue;
		}

		/* update bin count */
		nBinId = nPos/nBinSize;
		if( (nBinId < 0 ) || (nBinId >= vBinNum[nChrId]) )
		{
			printf("Error: CisChrome_GenomeProfile_CountBinReads_SingleFile, inconsistent bin number!\n");
			exit(EXIT_FAILURE);
		}

		/* '-' strand */
		if(nStrand == 1)
		{
			nk = nPos-nExtLen+1;
			nEndPos = nPos;
			nPos = nk;

			nEndBinId = nBinId;

			nBinId = nPos/nBinSize;
			if(nBinId < 0)
				nBinId = 0;
			if(nBinId >= vBinNum[nChrId])
				nBinId = vBinNum[nChrId]-1;
		}
		/* '+' strand */
		else
		{
			nEndPos = nPos+nExtLen-1;
			nEndBinId = nEndPos/nBinSize;
			if(nEndBinId < 0)
				nEndBinId = 0;
			if(nEndBinId >= vBinNum[nChrId])
				nEndBinId = vBinNum[nChrId]-1;
		}

		if(nEndPos < nPos)
		{
			printf("Error: CisChrome_GenomeProfile_CountBinReads_SingleFile, position order not correct!\n");
			exit(EXIT_FAILURE);
		}

		if(nBinId == nEndBinId)
		{
			fTemp = (float)(nEndPos-nPos+1)/(float)nBinSize;
			vBinC[nChrId][nBinId] += (float)(fTemp*dNormFactor);
		}
		else
		{
			fTemp = (float)(nBinSize - nPos%nBinSize)/(float)nBinSize;
			vBinC[nChrId][nBinId] += (float)(fTemp*dNormFactor);

			fTemp = (float)(nEndPos%nBinSize+1)/(float)nBinSize;
			vBinC[nChrId][nEndBinId] += (float)(fTemp*dNormFactor);

			for(nk=nBinId+1; nk<nEndBinId; nk++)
			{
				vBinC[nChrId][nk] += (float)(dNormFactor);
			}
		}
	}	

	/* close alignment file */
	fclose(fpRead);

	/* save & release memory */
	if(nExportBAR == 1)
	{
		sprintf(strOutFileName, "%s%s.bar.txt", strOutputPath, strOutputFile);
		fpRead = NULL;
		fpRead = fopen(strOutFileName, "w");
		if(fpRead == NULL)
		{
			printf("Error: cannot open output file!\n");
			exit(EXIT_FAILURE);
		}

		fprintf(fpRead, "#chr\tpos\t%s_binc\n", strFileName);
		fprintf(fpRead, "#chr\tpos\t1\n");

		/* process chromosome by chromosome */
		for(ni=0; ni<nChrNum; ni++)
		{
			nPos = nBinSize/2;
			for(nx=0; nx<vBinNum[ni]; nx++)
			{
				if(vBinC[ni][nx] > 1e-6)
					fprintf(fpRead, "%s\t%d\t%f\n", vChrName[ni]->m_pString, nPos, vBinC[ni][nx]);
				nPos += nBinSize;
			}
		}

		fclose(fpRead);
	}

	for(ni=0; ni<nChrNum; ni++)
	{
		sprintf(strOutFileName, "%s%s_%s.bincount", strOutputPath, strOutputFile, vChrName[ni]->m_pString);
		TileMapv2_SaveToBinaryFile((void *)(vBinC[ni]), sizeof(float), vBinNum[ni], strOutFileName);

		free(vBinC[ni]);
		vBinC[ni] = NULL;
	}
	free(vBinC);
	free(vBinNum);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_ComputeBinStats()                                            */
/*  Compute summary statistics for genomic bins.                           */
/* ----------------------------------------------------------------------- */ 
int CisChrome_ComputeBinStats(int nGroupNum, int nDatasetNum, int nSampleNum, 
		int nPaired, int nTwosided, int *vGroupId, int *vDatasetId, int **vRepNum,
		struct tagString **vSampleName, int nBinSize, int nWinSize,
		char strWorkFolder[], char strProjectTitle[], 
		int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen,
		int nGridNum, double *pMinT, double *pMaxT)
{
	/* define */
	int ni,nj,nk,nx,ny;
	int nBinNum;
	double dS2M = 0.0;
	double dS2S = 0.0;
	double dS2Scale = 0.0;
	double dTM = 0.0;
	double dTS = 0.0;
	int nTotalBinNum = 0;
	int nEffectBinNum = 0;
	int nDf;
	double dB = 0.0;
	double dMinT,dMaxT;
	double dStepsize;
	float *vThist;
	float fTCount;
	float *vTCDF;
	char strOutFileName[MED_LINE_LENGTH];
	
	/* process dataset by dataset to get initial statistics */
	for(nj=0; nj<nDatasetNum; nj++)
	{
		/* process chromosome by chromosome to get initial statistics */
		for(ni=0; ni<nChrNum; ni++)
		{
			nBinNum = pChrLen->pMatElement[ni]/nBinSize;
			if(pChrLen->pMatElement[ni] % nBinSize != 0)
				nBinNum++;

			CisChrome_ComputeBinStats_Chr_Initial(nGroupNum, nj,
				nSampleNum, nPaired, vGroupId, vDatasetId, vRepNum,
				vSampleName, nBinSize, nWinSize, nBinNum,
				strWorkFolder, strProjectTitle,
				vChrName[ni]->m_pString, 
				&dS2M, &dS2S, &nEffectBinNum);

			nTotalBinNum += nBinNum;
		}

		/* compute shrinkage factor */
		nDf = 0;
		if(nPaired == 0)
		{
			for(nk=0; nk<nGroupNum; nk++)
			{
				nDf += vRepNum[nk][nj]-1;
				if(vRepNum[nk][nj] > 0.5)
					dS2Scale += 1.0/vRepNum[nk][nj];
			}
		}
		else
		{
			nDf += vRepNum[0][nj]-1;
			if(vRepNum[0][nj] > 0.5)
				dS2Scale += 1.0/vRepNum[0][nj];
		}

		if(nDf > 0)
		{
			dS2M /= nEffectBinNum;
			dS2S -= nEffectBinNum*dS2M*dS2M;
			
			dB = 2.0*(nEffectBinNum-1)/(nDf+2.0)/nEffectBinNum + 2.0*dS2M*dS2M*(nEffectBinNum-1)/(nDf+2.0)/dS2S;
			if(dB < 0.0)
				dB = 0.0;
			if(dB > 1.0)
				dB = 1.0;
		}
		else
		{
			dS2M = 1.0;
			dS2S = 0.0;
			dS2Scale = 1.0;
			nDf = 0;
			dB = 1.0;
		}

		/* process chromosome by chromosome to get log2 fc and t-statistics with variance shrinkage */
		dMinT = 1e6;
		dMaxT = -1e6;
		for(ni=0; ni<nChrNum; ni++)
		{
			nBinNum = pChrLen->pMatElement[ni]/nBinSize;
			if(pChrLen->pMatElement[ni] % nBinSize != 0)
				nBinNum++;

			CisChrome_ComputeBinStats_Chr_Tstat(nj, strWorkFolder, strProjectTitle,  
				vChrName[ni]->m_pString, nBinNum, dB, dS2M, dS2Scale, &dMinT, &dMaxT);
		}

		/* get statistics distribution */
		vThist = NULL;
		vThist = (float *)calloc(nGridNum, sizeof(float));
		if(vThist == NULL)
		{
			printf("Error: CisChrome_ComputeBinStats, cannot allocate memory for the histogram!\n");
			exit(EXIT_FAILURE);
		}

		dStepsize = (dMaxT-dMinT)/nGridNum;
		pMinT[nj] = dMinT;
		pMaxT[nj] = dMaxT;

		for(ni=0; ni<nChrNum; ni++)
		{
			nBinNum = pChrLen->pMatElement[ni]/nBinSize;
			if(pChrLen->pMatElement[ni] % nBinSize != 0)
				nBinNum++;

			CisChrome_ComputeBinStats_Chr_Thist(nj, strWorkFolder, strProjectTitle,  
				vChrName[ni]->m_pString, nBinNum, nGridNum, dMinT, dMaxT, dStepsize, vThist);
		}

		vTCDF = NULL;
		vTCDF = (float *)calloc(nGridNum, sizeof(float));
		if(vTCDF == NULL)
		{
			printf("Error: CisChrome_ComputeBinStats, cannot allocate memory for the CDF histogram!\n");
			exit(EXIT_FAILURE);
		}

		fTCount = 0.0; 
		for(ni=0; ni<nGridNum; ni++)
		{
			fTCount += vThist[ni];
		}

		if(nTwosided == 0)
		{
			nx = nGridNum-1;
			vTCDF[nx] = vThist[nx];
			nx--;
			for(; nx>=0; nx--)
			{
				vTCDF[nx] = vTCDF[nx+1]+vThist[nx];
			}
		}
		else
		{
			nx = (int)((0-dMinT)/dStepsize);
			if(nx < 0)
				nx = 0;
			if(nx >= nGridNum)
				nx = nGridNum-1;

			ny = 0;
			vTCDF[ny] = vThist[ny];
			ny++;
			for(; ny<nx; ny++)
			{
				vTCDF[ny] = vTCDF[ny-1]+vThist[ny];
			}
			ny = nGridNum-1;
			vTCDF[ny] = vThist[ny];
			ny--;
			for(; ny>nx; ny--)
			{
				vTCDF[ny] = vTCDF[ny+1]+vThist[ny];
			}

			vTCDF[nx] = fTCount;
		}

		
		for(ni=0; ni<nGridNum; ni++)
		{
			vThist[ni] = (float)((vThist[ni]+1e-6)/(fTCount+nGridNum*1e-6));
			vTCDF[ni] = (float)((vTCDF[ni]+1e-6)/(fTCount+1e-6));
		}

		/* save data to files */
		sprintf(strOutFileName, "%s%s_%d_t.hist", strWorkFolder, strProjectTitle, nj);
		TileMapv2_SaveToBinaryFile((void *)vThist, sizeof(float), nGridNum, strOutFileName);

		sprintf(strOutFileName, "%s%s_%d_tcdf.hist", strWorkFolder, strProjectTitle, nj);
		TileMapv2_SaveToBinaryFile((void *)vTCDF, sizeof(float), nGridNum, strOutFileName);
		
		free(vThist);
		free(vTCDF);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_ComputeBinStats_Chr_Initial()                                */
/*  Compute initial summary statistics for genomic bins for a single       */
/*  chromosome.                                                            */
/* ----------------------------------------------------------------------- */ 
int CisChrome_ComputeBinStats_Chr_Initial(int nGroupNum, int nDatasetId,
			int nSampleNum, int nPaired, int *vGroupId, int *vDatasetId, 
			int **vRepNum, struct tagString **vSampleName, int nBinSize, 
			int nWinSize, int nBinNum, char strWorkFolder[], char strProjectTitle[],
			char strChr[], double *pS2M, double *pS2S, int *pEffectBinNum)
{
	/* define */
	int nHalfWinNum = 0;
	int nTemp;
	int nIPNum;
	int nCTNum;
	int nFileNum;
	int nGi,nDi;
	FILE **vfpIn;
	int ni,nj,nk,nl,nx;
	
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));

	/* data */
	float **vD = NULL;
	/* window average */
	float **vW = NULL;
	/* group mean */
	float *vM0 = NULL;
	float *vM1 = NULL;
	/* sample variance */
	float *vV = NULL;
	/* log2 fc */
	float *vFC = NULL;
	
	char strOutFileName[MED_LINE_LENGTH];
	char strInFileName[MED_LINE_LENGTH];
	int nWinNum;
	float dSum;
	float dTemp;
	int nDf;

	/* struct DOUBLEMATRIX *pLogNormFactor; */
	double dLog2 = log(2.0);
	double dLog10 = log(10.0);



	/* init */
	nHalfWinNum = (int)(nWinSize/nBinSize)-1;
	if(nHalfWinNum < 0)
		nHalfWinNum = 0;
	nTemp = nHalfWinNum%2;
	if(nTemp == 0) 
		nHalfWinNum = (int)(nHalfWinNum/2);
	else
		nHalfWinNum = (int)(nHalfWinNum/2)+1;
	nWinNum = 2*nHalfWinNum+1;
	
	if(nGroupNum == 1)
	{
		nIPNum = vRepNum[0][nDatasetId];
		nCTNum = 0;
	}
	else
	{
		nIPNum = vRepNum[0][nDatasetId];
		nCTNum = vRepNum[1][nDatasetId];

		if(nPaired == 1)
		{
			if(nIPNum != nCTNum)
			{
				printf("Error: CisChrome_ComputeBinStats_Chr_Initial, the number of samples in two conditions are not the same!\n");
				exit(EXIT_FAILURE);
			}
		}
	}

	nFileNum = nIPNum+nCTNum;

	vD = NULL;
	vD = (float **)calloc(nFileNum, sizeof(float *));
	if(vD == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Initial, cannot create vector for data!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nFileNum; ni++)
	{
		vD[ni] = (float *)calloc(nBinNum, sizeof(float));
		if(vD[ni] == NULL)
		{
			printf("Error: CisChrome_ComputeBinStats_Chr_Initial, cannot create memory block for data!\n");
			exit(EXIT_FAILURE);
		}
	}

	vW = NULL;
	vW = (float **)calloc(nFileNum, sizeof(float *));
	if(vW == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Initial, cannot create vector for window statistics!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nFileNum; ni++)
	{
		vW[ni] = (float *)calloc(nBinNum, sizeof(float));
		if(vW[ni] == NULL)
		{
			printf("Error: CisChrome_ComputeBinStats_Chr_Initial, cannot create memory block for data!\n");
			exit(EXIT_FAILURE);
		}
	}

	vM0 = (float *)calloc(nBinNum, sizeof(float));
	if(vM0 == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Initial, cannot allocate memory for M0!\n");
		exit(EXIT_FAILURE);
	}

	vM1 = (float *)calloc(nBinNum, sizeof(float));
	if(vM1 == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Initial, cannot allocate memory for M1!\n");
		exit(EXIT_FAILURE);
	}

	vV = (float *)calloc(nBinNum, sizeof(float));
	if(vV == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Initial, cannot allocate memory for V!\n");
		exit(EXIT_FAILURE);
	}

	vFC = (float *)calloc(nBinNum, sizeof(float));
	if(vFC == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Initial, cannot allocate memory for FC!\n");
		exit(EXIT_FAILURE);
	}

	vfpIn = NULL;
	vfpIn = (FILE **)calloc(nFileNum, sizeof(FILE *));
	if(vfpIn == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Initial, cannot create vector for file pointers!\n");
		exit(EXIT_FAILURE);
	}

	/* load data */
	nj = 0;
	nk = nIPNum;
	for(ni=0; ni<nSampleNum; ni++)
	{
		nGi = vGroupId[ni];
		nDi = vDatasetId[ni]-1;

		if(nDi != nDatasetId)
			continue;

		if(nGi == 1)
		{
			nl = nj;
			nj++;
		}
		else
		{
			nl = nk;
			nk++;
		}

		sprintf(strInFileName, "%s%s_%s.bincount", strWorkFolder, vSampleName[ni]->m_pString, strChr);
		vfpIn[nl] = fopen(strInFileName, "rb");
		if(vfpIn[nl] == NULL)
		{
			printf("Error: CisChrome_ComputeBinStats_Chr_Initial, cannot open file %s!\n", strInFileName);
			exit(EXIT_FAILURE);
		}
		
		if( little_endian_fread(vD[nl], sizeof(float), nBinNum, vfpIn[nl], little_endian_machine) != nBinNum)
		{
			printf("Error: CisChrome_ComputeBinStats_Chr_Initial, incorrect loading, number of bins inconsistent!\n");
			exit(EXIT_FAILURE);
		}

		fclose(vfpIn[nl]);
	}
	free(vfpIn);

	/* compute window average */
	for(ni=0; ni<nFileNum; ni++)
	{
		if(nBinNum <= nHalfWinNum)
			continue;
				
		dSum = 0.0;
		for(nj=0; nj<nHalfWinNum; nj++)
		{
			dSum += vD[ni][nj];
		}

		for(nj=0; nj<nBinNum; nj++)
		{
			nk = nj+nHalfWinNum;
			if(nk < nBinNum)
				dSum = dSum + vD[ni][nk];
			else
				nk = nBinNum-1;

			nl = nj-nHalfWinNum-1;
			if(nl >= 0)
				dSum = dSum - vD[ni][nl];
			else
				nl = -1;
			nx = nk-nl;

			vW[ni][nj] = (float)(log(1.0+dSum)/dLog2);

			if(nPaired == 0)
			{
				if(ni < nIPNum)
					vM1[nj] += vW[ni][nj];
				else
					vM0[nj] += vW[ni][nj];
			}
			else
			{
				if(ni >= nIPNum)
				{
					vW[ni-nIPNum][nj] -= vW[ni][nj];
					vM1[nj] += vW[ni-nIPNum][nj];
				}
			}
		}
	}

	/* compute mean and variance */
	nDf = 0;
	if(nPaired == 0)
	{
		if(nCTNum > 0)
			nDf += nCTNum-1;
		if(nIPNum > 0)
			nDf += nIPNum-1;
	}
	else
	{
		nDf = nIPNum-1;
	}

	for(nj=0; nj<nBinNum; nj++)
	{
		if(nPaired == 0)
		{
			if(nIPNum > 0)
				vM1[nj] /= nIPNum;

			if(nCTNum > 0)
				vM0[nj] /= nCTNum;

			vFC[nj] = vM1[nj]-vM0[nj];
			
		}
		else
		{
			if(nIPNum > 0)
				vM1[nj] /= nIPNum;
			vFC[nj] = vM1[nj];
		}

		if(nDf > 0)
		{
			for(ni=0; ni<nFileNum; ni++)
			{
				if(ni < nIPNum)
				{
					dTemp = vW[ni][nj] - vM1[nj];
					vV[nj] += dTemp*dTemp;
				}
				else
				{
					if(nPaired == 0)
					{
						dTemp = vW[ni][nj] - vM0[nj];
						vV[nj] += dTemp*dTemp;
					}
				}
			}

			vV[nj] /= nDf;

			if(vV[nj] > 1e-6)
			{
				*pS2M += vV[nj];
				*pS2S += (vV[nj]*vV[nj]);
				*pEffectBinNum = (*pEffectBinNum)+1;
			}
		}
	}

	/* save data to files */
	sprintf(strOutFileName, "%s%s_%d_%s.log2fc", strWorkFolder, strProjectTitle, nDatasetId, strChr);
	TileMapv2_SaveToBinaryFile((void *)vFC, sizeof(float), nBinNum, strOutFileName);

	sprintf(strOutFileName, "%s%s_%d_%s.v", strWorkFolder, strProjectTitle, nDatasetId, strChr);
	TileMapv2_SaveToBinaryFile((void *)vV, sizeof(float), nBinNum, strOutFileName);
	
	/* release memory */
	for(ni=0; ni<nFileNum; ni++)
	{
		free(vD[ni]);
		vD[ni] = NULL;

		free(vW[ni]);
		vW[ni] = NULL;
	}
	free(vD);
	free(vW);
	free(vM0);
	free(vM1);
	free(vV);
	free(vFC);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_ComputeBinStats_Chr_Tstat()                                  */
/*  Compute t-statistics.                                                  */
/* ----------------------------------------------------------------------- */ 
int CisChrome_ComputeBinStats_Chr_Tstat(int nDatasetId, char strWorkFolder[], 
			char strProjectTitle[], char strChr[], int nBinNum, double dB, 
			double dVM, double dVScale, double *pMinT, double *pMaxT)
{
	/* define */
	/* sample variance */
	float *vV = NULL;
	/* log2 fc */
	float *vM = NULL;
	/* t-stat */
	float *vT = NULL;
	/* file */
	FILE *fpIn;
	/* use var */
	int nUseVar = 0;

	char strOutFileName[MED_LINE_LENGTH];
	char strInFileName[MED_LINE_LENGTH];
	double dTemp;

	int ni;
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));

	/* init */
	vM = (float *)calloc(nBinNum, sizeof(float));
	if(vM == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Tstatt, cannot allocate memory for M!\n");
		exit(EXIT_FAILURE);
	}

	vV = (float *)calloc(nBinNum, sizeof(float));
	if(vV == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Tstat, cannot allocate memory for V!\n");
		exit(EXIT_FAILURE);
	}

	vT = (float *)calloc(nBinNum, sizeof(float));
	if(vT == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Tstat, cannot allocate memory for T!\n");
		exit(EXIT_FAILURE);
	}

	/* load D */
	sprintf(strInFileName, "%s%s_%d_%s.log2fc", strWorkFolder, strProjectTitle, nDatasetId, strChr);
	fpIn = NULL;
	fpIn = fopen(strInFileName, "rb");
	if(fpIn == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Tstat, cannot open file %s!\n", strInFileName);
		exit(EXIT_FAILURE);
	}
	
	if( little_endian_fread(vM, sizeof(float), nBinNum, fpIn, little_endian_machine) != nBinNum)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Tstat, incorrect loading, number of bins inconsistent!\n");
		exit(EXIT_FAILURE);
	}

	fclose(fpIn);

	/* load V */
	sprintf(strInFileName, "%s%s_%d_%s.v", strWorkFolder, strProjectTitle, nDatasetId, strChr);
	fpIn = NULL;
	fpIn = fopen(strInFileName, "rb");
	if(fpIn == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Tstat, cannot open file %s!\n", strInFileName);
		exit(EXIT_FAILURE);
	}
	
	if( little_endian_fread(vV, sizeof(float), nBinNum, fpIn, little_endian_machine) != nBinNum)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Tstat, incorrect loading, number of bins inconsistent!\n");
		exit(EXIT_FAILURE);
	}

	fclose(fpIn);

	/* compute t-stat */
	for(ni=0; ni<nBinNum; ni++)
	{
		dTemp = (1.0-dB)*vV[ni]+dB*dVM;
		dTemp = sqrt(dVScale*dTemp)+1e-6;
		vT[ni] = (float)(vM[ni]/dTemp);

		if(vT[ni] < (*pMinT))
			(*pMinT) = vT[ni];
		if(vT[ni] > (*pMaxT))
			(*pMaxT) = vT[ni];
	}

	/* save t-stat */
	sprintf(strOutFileName, "%s%s_%d_%s.t", strWorkFolder, strProjectTitle, nDatasetId, strChr);
	TileMapv2_SaveToBinaryFile((void *)vT, sizeof(float), nBinNum, strOutFileName);

	/* release memory */
	free(vM);
	free(vV);
	free(vT);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_ComputeBinStats_Chr_Thist()                                  */
/*  Get distribution of t-statistics.                                      */
/* ----------------------------------------------------------------------- */ 
int CisChrome_ComputeBinStats_Chr_Thist(int nDatasetId, char strWorkFolder[], 
		char strProjectTitle[],  char strChr[], int nBinNum, int nGridNum,
		double dMinT, double dMaxT, double dStep, float *vThist)
{
	/* define */
	/* sample variance */
	float *vV = NULL;
	/* log2 fc */
	float *vM = NULL;
	/* t-stat */
	float *vT = NULL;
	/* file */
	FILE *fpIn;

	char strInFileName[MED_LINE_LENGTH];

	int ni,nId;
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));

	/* init */
	vM = (float *)calloc(nBinNum, sizeof(float));
	if(vM == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Thist, cannot allocate memory for M!\n");
		exit(EXIT_FAILURE);
	}

	vV = (float *)calloc(nBinNum, sizeof(float));
	if(vV == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Thist, cannot allocate memory for V!\n");
		exit(EXIT_FAILURE);
	}

	vT = (float *)calloc(nBinNum, sizeof(float));
	if(vT == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Thist, cannot allocate memory for T!\n");
		exit(EXIT_FAILURE);
	}

	/* load D */
	sprintf(strInFileName, "%s%s_%d_%s.log2fc", strWorkFolder, strProjectTitle, nDatasetId, strChr);
	fpIn = NULL;
	fpIn = fopen(strInFileName, "rb");
	if(fpIn == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Thist, cannot open file %s!\n", strInFileName);
		exit(EXIT_FAILURE);
	}
	
	if( little_endian_fread(vM, sizeof(float), nBinNum, fpIn, little_endian_machine) != nBinNum)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Thist, incorrect loading, number of bins inconsistent!\n");
		exit(EXIT_FAILURE);
	}

	fclose(fpIn);

	/* load V */
	sprintf(strInFileName, "%s%s_%d_%s.v", strWorkFolder, strProjectTitle, nDatasetId, strChr);
	fpIn = NULL;
	fpIn = fopen(strInFileName, "rb");
	if(fpIn == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Thist, cannot open file %s!\n", strInFileName);
		exit(EXIT_FAILURE);
	}
	
	if( little_endian_fread(vV, sizeof(float), nBinNum, fpIn, little_endian_machine) != nBinNum)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Thist, incorrect loading, number of bins inconsistent!\n");
		exit(EXIT_FAILURE);
	}

	fclose(fpIn);

	/* load T */
	sprintf(strInFileName, "%s%s_%d_%s.t", strWorkFolder, strProjectTitle, nDatasetId, strChr);
	fpIn = NULL;
	fpIn = fopen(strInFileName, "rb");
	if(fpIn == NULL)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Thist, cannot open file %s!\n", strInFileName);
		exit(EXIT_FAILURE);
	}
	
	if( little_endian_fread(vT, sizeof(float), nBinNum, fpIn, little_endian_machine) != nBinNum)
	{
		printf("Error: CisChrome_ComputeBinStats_Chr_Thist, incorrect loading, number of bins inconsistent!\n");
		exit(EXIT_FAILURE);
	}

	fclose(fpIn);

	/* compute t-stat */
	for(ni=0; ni<nBinNum; ni++)
	{
		if( (fabs(vM[ni])<1e-6) && (fabs(vV[ni])<1e-6) )
			continue;


		nId = (int)((vT[ni]-dMinT)/dStep);
		if(nId >= nGridNum)
			nId = nGridNum-1;

		vThist[nId] += 1.0;
	}

	/* release memory */
	free(vM);
	free(vV);
	free(vT);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_FindActiveTFBS()                                             */
/*  Find active transcription factor binding sites.                        */
/* ----------------------------------------------------------------------- */ 
int CisChrome_FindActiveTFBS(char strWorkFolder[], char strProjectTitle[],
		int nDatasetNum, int nBinSize, int nWinSize, int nTwoSided,
		int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen, 
		int nGridNum, char strMotifFolder[], int nMotifNum, 
		struct tagString **vTFName, struct tagString **vMotifFilePath,
		float **vMotifStat, float ** vMotifCDF, double *vMinT, double *vMaxT, 
		double dFDRCut, int *vTotalSiteNum, int *vActiveSiteNum, 
		double *vActiveSitePrc, double *vTFExpress, 
		struct tagRefGene **vTFRefGene, int nUp, int nDown, int nTSS)
{
	/* define */
	int ni,nj,nk;
	char strMotifFileName[MED_LINE_LENGTH];
	int nTFBSNum;
	struct tagString **vTFBSChr;
	int *vTFBSPos1;
	int *vTFBSPos2;
	unsigned char *vTFBSStrand;
	double dStepSize;
	int nActiveTFBSNum;
	int nTFPos1,nTFPos2;

	/* process motif by motif */
	printf("\n");
	printf("/* ------------------------------------------------- */\n");
	printf("/* Find active TFBS                                  */\n");
	printf("/* ------------------------------------------------- */\n");
		
	for(ni=0; ni<nMotifNum; ni++)
	{
		/* init */
		nTFBSNum = 0;
		vTFBSChr = NULL;
		vTFBSPos1 = NULL;
		vTFBSPos2 = NULL;
		vTFBSStrand = NULL;
		
		if(nTSS == 1)
		{
			if(vTFRefGene[ni]->chStrand == '-')
			{
				nTFPos1 = vTFRefGene[ni]->nTxEnd-nDown;
				nTFPos2 = vTFRefGene[ni]->nTxEnd+nUp;
			}
			else
			{
				nTFPos1 = vTFRefGene[ni]->nTxStart-nUp;
				nTFPos2 = vTFRefGene[ni]->nTxStart+nDown;
			}
		}
		else
		{
			if(vTFRefGene[ni]->chStrand == '-')
			{
				nTFPos1 = vTFRefGene[ni]->nTxStart-nDown;
				nTFPos2 = vTFRefGene[ni]->nTxEnd+nUp;
			}
			else
			{
				nTFPos1 = vTFRefGene[ni]->nTxStart-nUp;
				nTFPos2 = vTFRefGene[ni]->nTxEnd+nDown;
			}
		}

		/* get motif site number */
		sprintf(strMotifFileName, "%s%s", strMotifFolder, vMotifFilePath[ni]->m_pString);
		printf("%d: Processing %s ...",  ni+1, vMotifFilePath[ni]->m_pString);

		/* read motif sites */
		CisChrome_FindActiveTFBS_LoadTFBS(strMotifFileName, &nTFBSNum, &vTFBSChr, 
			&vTFBSPos1, &vTFBSPos2, &vTFBSStrand);
		vTotalSiteNum[ni] = nTFBSNum;

		/* process dataset by dataset */
		for(nj=0; nj<nDatasetNum; nj++)
		{
			nk = ni*nDatasetNum+nj;
			dStepSize = (vMaxT[nj]-vMinT[nj])/nGridNum;
			nActiveTFBSNum = CisChrome_FindActiveTFBS_Call(strWorkFolder, strProjectTitle, nj, 
				(ni+1), vTFName[ni]->m_pString, nTFBSNum, vTFBSChr, vTFBSPos1, vTFBSPos2, 
				vTFBSStrand, nChrNum, vChrName, pChrLen, nBinSize, nWinSize, nTwoSided,
				(vMotifStat[ni]+nj*nGridNum), (vMotifCDF[ni]+nj*nGridNum), vMinT[nj], vMaxT[nj], 
				dStepSize, nGridNum, dFDRCut, vTFRefGene[ni]->strChrom, nTFPos1, nTFPos2, 
				(vTFExpress+nk));
						
			vActiveSiteNum[nk] = nActiveTFBSNum;
			vActiveSitePrc[nk] = (nActiveTFBSNum+1e-6)/(nTFBSNum+1e-6);
		}

		/* release memory */
		for(nj=0; nj<nTFBSNum; nj++)
		{
			DeleteString(vTFBSChr[nj]);
			vTFBSChr[nj] = NULL;
		}
		free(vTFBSPos1);
		free(vTFBSPos2);
		free(vTFBSStrand);
		free(vTFBSChr);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_FindActiveTFBS_LoadTFBS()                                    */
/*  Load motif sites.                                                      */
/* ----------------------------------------------------------------------- */ 
int CisChrome_FindActiveTFBS_LoadTFBS(char strMotifFileName[], int *pTFBSNum, 
		struct tagString ***pvTFBSChr, int **pvTFBSPos1, int **pvTFBSPos2,
		unsigned char **pvTFBSStrand)
{
	/* define */
	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	char strAlias[MED_LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nPos1,nPos2;
	char chStrand;
	int ni;

	/* init */
	*pTFBSNum = 0;
	fpIn = NULL;
	fpIn = fopen(strMotifFileName, "rt");
	if(fpIn == NULL)
	{
		printf("Error: CisChrome_FindActiveTFBS_LoadTFBS, cannot open motif site file!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;
		
		*pTFBSNum += 1;
	}
	fclose(fpIn);
	printf("motif_site_num=%d\n", *pTFBSNum);

	/* create memory */
	*pvTFBSChr = NULL;
	*pvTFBSChr = (struct tagString **)calloc(*pTFBSNum, sizeof(struct tagString *));
	if(*pvTFBSChr == NULL)
	{
		printf("Error: CisChrome_FindActiveTFBS_LoadTFBS, cannot allocate memory for motif site chromosomes!\n");
		exit(EXIT_FAILURE);
	}

	*pvTFBSPos1 = NULL;
	*pvTFBSPos1 = (int *)calloc(*pTFBSNum, sizeof(int));
	if(*pvTFBSPos1 == NULL)
	{
		printf("Error: CisChrome_FindActiveTFBS_LoadTFBS, cannot allocate memory for motif site locations!\n");
		exit(EXIT_FAILURE);
	}

	*pvTFBSPos2 = NULL;
	*pvTFBSPos2 = (int *)calloc(*pTFBSNum, sizeof(int));
	if(*pvTFBSPos2 == NULL)
	{
		printf("Error: CisChrome_FindActiveTFBS_LoadTFBS, cannot allocate memory for motif site locations!\n");
		exit(EXIT_FAILURE);
	}

	*pvTFBSStrand = NULL;
	*pvTFBSStrand = (unsigned char *)calloc(*pTFBSNum, sizeof(unsigned char));
	if(*pvTFBSStrand == NULL)
	{
		printf("Error: CisChrome_FindActiveTFBS_LoadTFBS, cannot allocate memory for motif site strands!\n");
		exit(EXIT_FAILURE);
	}

	/* read motif sites */
	fpIn = NULL;
	fpIn = fopen(strMotifFileName, "rt");
	if(fpIn == NULL)
	{
		printf("Error: CisChrome_FindActiveTFBS_LoadTFBS, cannot open motif site file!\n");
		exit(EXIT_FAILURE);
	}
	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %s %d %d %c", strAlias, strChr, &nPos1, &nPos2, &chStrand);

		StringAddTail((*pvTFBSChr)+ni, strChr);
		(*pvTFBSPos1)[ni] = nPos1;
		(*pvTFBSPos2)[ni] = nPos2;
		
		if( (chStrand == '-') || (chStrand == 'R') ) 
			(*pvTFBSStrand)[ni] = 1;
		else
			(*pvTFBSStrand)[ni] = 0;

		ni++;
	}
	fclose(fpIn);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_FindActiveTFBS_Call()                                        */
/*  Call active motif sites.                                               */
/* ----------------------------------------------------------------------- */ 
int CisChrome_FindActiveTFBS_Call(char strWorkFolder[], char strProjectTitle[], 
		int nDatasetId, int nTFId, char strTFName[], int nTFBSNum, 
		struct tagString **vTFBSChr, int *vTFBSPos1, int *vTFBSPos2, 
		unsigned char *vTFBSStrand,	int nChrNum, struct tagString **vChrName, 
		struct INTMATRIX *pChrLen, int nBinSize, int nWinSize, int nTwoSided, 
		float *vMotifHist, float *vMotifCDF, double dMinT, double dMaxT, 
		double dStepSize, int nGridNum, double dFDRCut,
		char strTFChr[], int nTFPos1, int nTFPos2, double *pTFExpress)
{
	/* define */
	float **vData;
	int ni,nj,nChrId,nk,nPos;
	int nBinNum;
	char strInFileName[MED_LINE_LENGTH];
	char strOutFileName[MED_LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	struct DOUBLEMATRIX *pTFBSScore;
	int nActiveTFBSNum = 0;
	double dSum;
	float *vBGHist;
	float *vBGCDF;
	float fCutL,fCutH;
	int nZId;
	double dTemp;

	/* load data */
	vData = NULL;
	vData = (float **)calloc(nChrNum, sizeof(float *));
	if(vData == NULL)
	{
		printf("Error: CisChrome_FindActiveTFBS_Call, cannot allocate memory for the enrichment data!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nChrNum; ni++)
	{
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		vData[ni] = (float *)calloc(nBinNum, sizeof(float));
		if(vData[ni] == NULL)
		{
			printf("Error: CisChrome_FindActiveTFBS_Call, cannot allocate memory for the enrichment data!\n");
			exit(EXIT_FAILURE);
		}

		sprintf(strInFileName, "%s%s_%d_%s.t", strWorkFolder, strProjectTitle, nDatasetId, vChrName[ni]->m_pString);
		fpIn = NULL;
		fpIn = fopen(strInFileName, "rb");
		if(fpIn == NULL)
		{
			printf("Error: CisChrome_FindActiveTFBS_Call, cannot open file %s!\n", strInFileName);
			exit(EXIT_FAILURE);
		}
	
		if( little_endian_fread(vData[ni], sizeof(float), nBinNum, fpIn, little_endian_machine) != nBinNum)
		{
			printf("Error: CisChrome_FindActiveTFBS_Call, incorrect loading, number of bins inconsistent!\n");
			exit(EXIT_FAILURE);
		}

		fclose(fpIn);
	}

	/* get enrichment score for the TF */
	nChrId = SeqClust_FindChr(vChrName, nChrNum, strTFChr);
	if( (nChrId < 0) || (nChrId >= nChrNum) )
	{
		printf("Warning: CisChrome_FindActiveTFBS_Call, cannot find the matching chromosome for TF %d:%d!\n",
			nTFId, strTFName);
	}
	else if( (nTFPos2 <= 0) || (nTFPos1 >= pChrLen->pMatElement[nChrId]) )
	{
		printf("Warning: CisChrome_FindActiveTFBS_Call, the coordinates of TF %d:%d are beyond the chromosome length!\n",
			nTFId, strTFName);
	}
	else
	{
		if(nTFPos1 < 0)
			nTFPos1 = 0;
		if(nTFPos2 >= pChrLen->pMatElement[nChrId])
			nTFPos2 = pChrLen->pMatElement[nChrId]-1;

		ni = nTFPos1/nBinSize;
		nj = nTFPos2/nBinSize;

		dSum = 0.0;
		for(nk=ni; nk<=nj; nk++)
		{
			dSum += vData[nChrId][nk];
		}
		*pTFExpress = dSum/(nj-ni+1);
	}

	/* get enrichment score at each motif site */
	pTFBSScore = NULL;
	pTFBSScore = CreateDoubleMatrix(1, nTFBSNum);
	if(pTFBSScore == NULL)
	{
		printf("Error: CisChrome_FindActiveTFBS_Call, cannot create memory for TFBS scores!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nTFBSNum; ni++)
	{
		/* get nChrId */
		nChrId = SeqClust_FindChr(vChrName, nChrNum, vTFBSChr[ni]->m_pString);
		nPos = (vTFBSPos1[ni]+vTFBSPos2[ni])/2;

		if( (nChrId < 0) || (nChrId >= nChrNum) )
		{
			printf("Warning: CisChrome_FindActiveTFBS_Call, cannot find the matching chromosome for %s:%d-%d!\n",
				vTFBSChr[ni]->m_pString, vTFBSPos1[ni], vTFBSPos2[ni]);
			continue;
		}

		if(nPos >= pChrLen->pMatElement[nChrId])
		{
			/* printf("Warning: CisChrome_FindActiveTFBS_Call, motif site %s:%d-%d beyond the range!\n", 
				vTFBSChr[ni]->m_pString, vTFBSPos1[ni], vTFBSPos2[ni]); */
			continue;
		}

		nj = nPos/nBinSize;
		pTFBSScore->pMatElement[ni] = vData[nChrId][nj];
		nk = (int)((pTFBSScore->pMatElement[ni]-dMinT)/dStepSize);
		if(nk >= nGridNum)
			nk = nGridNum-1;
		vMotifHist[nk] += 1.0;
	}

	/* compute FDR */
	vBGHist = NULL;
	vBGHist = (float *)calloc(nGridNum, sizeof(float));
	if(vBGHist == NULL)
	{
		printf("Error: CisChrome_FindActiveTFBS_Call, cannot create memory for background T histogram!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strInFileName, "%s%s_%d_t.hist", strWorkFolder, strProjectTitle, nDatasetId);
	TileMapv2_LoadFromBinaryFile(vBGHist, sizeof(float), nGridNum, strInFileName);

	vBGCDF = NULL;
	vBGCDF = (float *)calloc(nGridNum, sizeof(float));
	if(vBGCDF == NULL)
	{
		printf("Error: CisChrome_FindActiveTFBS_Call, cannot create memory for background T CDF histogram!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strInFileName, "%s%s_%d_tcdf.hist", strWorkFolder, strProjectTitle, nDatasetId);
	TileMapv2_LoadFromBinaryFile(vBGCDF, sizeof(float), nGridNum, strInFileName);

	dSum = 0.0;
	for(ni=0; ni<nGridNum; ni++)
	{
		dSum += vMotifHist[ni];
	}

	if(nTwoSided == 0)
	{
		ni = nGridNum-1;
		vMotifCDF[ni] = vMotifHist[ni];
		ni--;
		for(; ni>=0; ni--)
		{
			vMotifCDF[ni] = vMotifCDF[ni+1]+vMotifHist[ni];
		}

		vMotifCDF[0] = (float)((vMotifCDF[0]+1e-6)/(dSum+1e-6)/vBGCDF[0]);
		for(ni=1; ni<nGridNum; ni++)
		{
			vMotifCDF[ni] = (float)((vMotifCDF[ni]+1e-6)/(dSum+1e-6)/vBGCDF[ni]);
			/* if(vMotifCDF[ni] < vMotifCDF[ni-1])
				vMotifCDF[ni] = vMotifCDF[ni-1]; */
		}
	}
	else
	{
		nZId = (int)((0.0-dMinT)/dStepSize);
		if(nZId < 0)
			nZId = 0;
		if(nZId >= nGridNum)
			nZId = nGridNum-1;

		ni = 0;
		vMotifCDF[ni] = vMotifHist[ni];
		ni++;
		for(; ni<nZId; ni++)
		{
			vMotifCDF[ni] = vMotifCDF[ni-1]+vMotifHist[ni];
		}
		
		ni = nGridNum-1;
		vMotifCDF[ni] = vMotifHist[ni];
		ni--;
		for(; ni>nZId; ni--)
		{
			vMotifCDF[ni] = vMotifCDF[ni+1]+vMotifHist[ni];
		}
		vMotifCDF[nZId] = (float)dSum;

		vMotifCDF[nZId] = (float)((vMotifCDF[nZId]+1e-6)/(dSum+1e-6)/vBGCDF[nZId]);
		for(ni=nZId+1; ni<nGridNum; ni++)
		{
			vMotifCDF[ni] = (float)((vMotifCDF[ni]+1e-6)/(dSum+1e-6)/vBGCDF[ni]);
			/* if(vMotifCDF[ni] < vMotifCDF[ni-1])
				vMotifCDF[ni] = vMotifCDF[ni-1]; */
		}
		for(ni=nZId-1; ni>=0; ni--)
		{
			vMotifCDF[ni] = (float)((vMotifCDF[ni]+1e-6)/(dSum+1e-6)/vBGCDF[ni]);
			/* if(vMotifCDF[ni] < vMotifCDF[ni+1])
				vMotifCDF[ni] = vMotifCDF[ni+1]; */
		}
	}

	for(ni=0; ni<nGridNum; ni++)
	{
		vMotifHist[ni] = (float)((vMotifHist[ni]+1e-6)/(dSum+nGridNum*1e-6));
		vMotifHist[ni] = (float)(vMotifHist[ni]/vBGHist[ni]);
	}

	fCutH = (float)dMaxT;
	fCutL = (float)dMinT;
	if(nTwoSided == 0)
	{
		for(ni=0; ni<nGridNum; ni++)
		{
			dTemp = 1.0/vMotifCDF[ni];
			if(dTemp <= dFDRCut)
			{
				fCutH = (float)(dMinT+dStepSize*ni);
				break;
			}
		}
	}
	else
	{
		nZId = (int)((0.0-dMinT)/dStepSize);
		if(nZId < 0)
			nZId = 0;
		if(nZId >= nGridNum)
			nZId = nGridNum-1;

	
		for(ni=nZId; ni<nGridNum; ni++)
		{
			dTemp = 1.0/vMotifCDF[ni];
			if(dTemp <= dFDRCut)
			{
				fCutH = (float)(dMinT+dStepSize*ni);
				break;
			}
		}
		for(ni=nZId; ni>=0; ni--)
		{
			dTemp = 1.0/vMotifCDF[ni];
			if(dTemp <= dFDRCut)
			{
				fCutL = (float)(dMinT+dStepSize*(ni+1));
				break;
			}
		}
	}

	/* export top motif sites */
	sprintf(strOutFileName, "%s%s_%d_%s_%d_ActiveTFBS.txt", strWorkFolder, strProjectTitle, nTFId, strTFName, nDatasetId);
	fpOut = NULL;
	fpOut = fopen(strOutFileName, "wt");
	if(fpOut == NULL)
	{
		printf("Error: CisChrome_FindActiveTFBS_Call, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fpOut, "#motifsite_id\tchr\tstart\tend\tstrand\tscore\n");
	for(ni=0; ni<nTFBSNum; ni++)
	{
		if(pTFBSScore->pMatElement[ni] >= fCutH)
		{
			fprintf(fpOut, "%d\t%s\t%d\t%d\t", ni, vTFBSChr[ni]->m_pString, vTFBSPos1[ni], vTFBSPos2[ni]);
			if(vTFBSStrand[ni] == 0)
				fprintf(fpOut, "+");
			else
				fprintf(fpOut, "-");

			fprintf(fpOut, "\t%f\n", pTFBSScore->pMatElement[ni]);
			nActiveTFBSNum++;
		}

		if(nTwoSided == 1)
		{
			if(pTFBSScore->pMatElement[ni] <= fCutL)
			{
				fprintf(fpOut, "%d\t%s\t%d\t%d\t", ni, vTFBSChr[ni]->m_pString, vTFBSPos1[ni], vTFBSPos2[ni]);
				if(vTFBSStrand[ni] == 0)
					fprintf(fpOut, "+");
				else
					fprintf(fpOut, "-");

				fprintf(fpOut, "\t%f\n", pTFBSScore->pMatElement[ni]);
				nActiveTFBSNum++;
			}
		}
	}
	fclose(fpOut);

	/* release memory */
	for(ni=0; ni<nChrNum; ni++)
	{
		free(vData[ni]);
		vData[ni] = NULL;
	}
	free(vData);
	DestroyDoubleMatrix(pTFBSScore);
	free(vBGHist);
	free(vBGCDF);

	/* return */
	return nActiveTFBSNum;
}

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_Export()                                                     */
/*  Export results.                                                        */
/* ----------------------------------------------------------------------- */ 
int CisChrome_Export(char strWorkFolder[], char strProjectTitle[],
		int nMotifNum, int nDatasetNum, int nGridNum, struct tagString **vTFName, 
		int *vTotalSiteNum, int *vActiveSiteNum, double *vActiveSitePrc,
		float **vMotifStat, double *vMinT, double *vMaxT, 
		double *vTFExpress, int nOutputType)
{
	/* define */
	FILE *fpOut;
	char strFileName[MED_LINE_LENGTH];
	int ni,nj,nk,nl;
	double dStep;
	double dLog2 = log(2.0);

	/* read motif sites */
	if(nOutputType == 0)
		sprintf(strFileName, "%s%s_dSummary.txt", strWorkFolder, strProjectTitle);
	else
		sprintf(strFileName, "%s%s_CDFSummary.txt", strWorkFolder, strProjectTitle);

	fpOut = NULL;
	fpOut = fopen(strFileName, "wt");
	if(fpOut == NULL)
	{
		printf("Error: CisChrome_Export, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#MotifID\tTFName\tTotal_Site_Num");
	for(nj=0; nj<nDatasetNum; nj++)
	{
		fprintf(fpOut, "\tD%d_TFExpr\tD%d_Active_Site_Num\tD%d_Active_Site_Prc", nj+1, nj+1, nj+1);
		for(nl=0; nl<nGridNum; nl++)
		{
			fprintf(fpOut, "\tD%d_Grid%d", nj+1, nl+1);
		}
	}
	fprintf(fpOut, "\n");

	fprintf(fpOut, "#MotifID\tTFName\tTotal_Site_Num");
	for(nj=0; nj<nDatasetNum; nj++)
	{
		dStep = (vMaxT[nj]-vMinT[nj])/nGridNum;
		
		fprintf(fpOut, "\tD%d_TFExpr\tD%d_Active_Site_Num\tD%d_Active_Site_Prc", nj+1, nj+1, nj+1);
		for(nl=0; nl<nGridNum; nl++)
		{
			fprintf(fpOut, "\t[%f,%f]", vMinT[nj]+dStep*nl, vMinT[nj]+dStep*(nl+1));
		}
	}
	fprintf(fpOut, "\n");

	fprintf(fpOut, "#MotifID\tTFName\tTotal_Site_Num");
	for(nj=0; nj<nDatasetNum; nj++)
	{
		dStep = (vMaxT[nj]-vMinT[nj])/nGridNum;
		
		fprintf(fpOut, "\tD%d_TFExpr\tD%d_Active_Site_Num\tD%d_Active_Site_Prc", nj+1, nj+1, nj+1);
		for(nl=0; nl<nGridNum; nl++)
		{
			fprintf(fpOut, "\t%f", vMinT[nj]+dStep*nl+dStep/2);
		}
	}
	fprintf(fpOut, "\n");
	
	for(ni=0; ni<nMotifNum; ni++)
	{
		fprintf(fpOut, "%d\t%s\t%d", ni+1, vTFName[ni]->m_pString, vTotalSiteNum[ni]);
		for(nj=0; nj<nDatasetNum; nj++)
		{
			nk = ni*nDatasetNum+nj;
			fprintf(fpOut, "\t%f\t%d\t%f", vTFExpress[nk], vActiveSiteNum[nk], vActiveSitePrc[nk]);
			
			nk = (nj+1)*nGridNum;
			for(nl=nj*nGridNum; nl<nk; nl++)
			{
				fprintf(fpOut, "\t%f", log(vMotifStat[ni][nl])/dLog2);
			}
		}
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}
