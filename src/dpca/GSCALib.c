/* ----------------------------------------------------------------------- */
/*  AffyLib.c : implementation of the affymetrix library                   */
/*  Author : Ji HongKai ; Time: 2004.08                                    */
/* ----------------------------------------------------------------------- */

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "limits.h"

#include "MathLib.h"
#include "MatrixLib.h"
#include "StringLib.h"
#include "GSCALib.h"
#include "GenomeLib.h"
#include "AffyLib.h"
#include "TilingArrayLib.h"

/* ----------------------------------------------------------------------- */ 
/*                            Global Variables                             */
/* ----------------------------------------------------------------------- */ 

/* ----------------------------------------------------------------------- */ 
/*                             Functions                                   */
/* ----------------------------------------------------------------------- */ 

/* ----------------------------------------------------------------------- */ 
/*  TexGenome_BuildIndex()                                                 */
/* ----------------------------------------------------------------------- */ 
int TexGenome_BuildIndex(char strSampleIndexFile[], char strDocIndexFile[], 
						 char strWordMapFile[], char strExportFolder[])
{
	/*char strSampleIndexFile[MED_LINE_LENGTH];
	char strDocIndexFile[MED_LINE_LENGTH];
	char strWordMapFile[MED_LINE_LENGTH];
	char strExportFolder[MED_LINE_LENGTH]; 
	*/
	
	/* define */
	int nSampleNum = 0;
	int nDocNum = 0;
	int nWordNum = 0;
	
	struct tagString **vDocIndex;
	float *vDocScore;

	struct tagString **vGSMIndex;
	float *vSampleScore;

	int *vSampleDocMap;


	char *vString;
	FILE *fpIn;
	int nLineLength;
	char *chp1,*chp2;
	char separators[]="\t ";
	char strWord[LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	int ni,nj,nk,nl;

	/* init */
	/* strcpy(strSampleIndexFile, "H:\\Projects\\JiHK_Lab\\NLP\\GPL1261_sample_index.txt");
	strcpy(strDocIndexFile, "H:\\Projects\\JiHK_Lab\\NLP\\all_doc_index.txt");
	strcpy(strWordMapFile, "H:\\Projects\\JiHK_Lab\\NLP\\tfidf_matrix.txt");
	strcpy(strExportFolder, "H:\\Projects\\JiHK_Lab\\NLP\\tfidf\\"); */

	AdjustDirectoryPath(strExportFolder);
	
	/* ----- STEP 1 -------- */
	/* create output folders */
	for(ni='a'; ni<='z'; ni++)
	{
		sprintf(strLine, "mkdir %s%c", strExportFolder, ni);
		system(strLine);
	}
	for(ni='0'; ni<='9'; ni++)
	{
		sprintf(strLine, "mkdir %s%c", strExportFolder, ni);
		system(strLine);
	}
	sprintf(strLine, "mkdir %sother", strExportFolder);
	system(strLine);

	/* ----- STEP 2 -------- */
	/* read sample index */
	fpIn = NULL;
	fpIn = fopen(strSampleIndexFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: TexGenome_BuildIndex, cannot open sample index file!\n");
		exit(EXIT_FAILURE);
	}

	nSampleNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nSampleNum++;
	}
	fclose(fpIn);

	if(nSampleNum <= 0)
	{
		printf("Error: TexGenome_BuildIndex, no sample available!");
		return PROC_SUCCESS;
	}

	vGSMIndex = NULL;
	vGSMIndex = (struct tagString **)calloc(nSampleNum, sizeof(struct tagString *));
	if(vGSMIndex == NULL)
	{
		printf("Error: TexGenome_BuildIndex, cannot create GSM index!");
		exit(EXIT_FAILURE);
	}


	fpIn = NULL;
	fpIn = fopen(strSampleIndexFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: TexGenome_BuildIndex, cannot open sample index file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		chp1 = strchr(strLine, '\t');
		if(chp1 != NULL)
			*chp1 = '\0';
		StringAddTail(vGSMIndex+ni, strLine);

		ni++;
	}
	fclose(fpIn);

	if(ni != nSampleNum)
	{
		printf("Error: TexGenome_BuildIndex, inconsistent sample number!\n");
		exit(EXIT_FAILURE);
	}

	/* ----- STEP 3 -------- */
	/* read document index */
	fpIn = NULL;
	fpIn = fopen(strDocIndexFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: TexGenome_BuildIndex, cannot open document index file!\n");
		exit(EXIT_FAILURE);
	}

	nDocNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nDocNum++;
	}
	fclose(fpIn);

	if(nDocNum <= 0)
	{
		printf("Error: TexGenome_BuildIndex, no document available!");
		return PROC_SUCCESS;
	}

	vDocIndex = NULL;
	vDocIndex = (struct tagString **)calloc(nDocNum, sizeof(struct tagString *));
	if(vDocIndex == NULL)
	{
		printf("Error: TexGenome_BuildIndex, cannot create document index!");
		exit(EXIT_FAILURE);
	}


	fpIn = NULL;
	fpIn = fopen(strDocIndexFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: TexGenome_BuildIndex, cannot open document index file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		chp1 = strchr(strLine, ':');
		if(chp1 != NULL)
			*chp1 = '\0';
		StringAddTail(vDocIndex+ni, strLine);

		ni++;
	}
	fclose(fpIn);

	if(ni != nDocNum)
	{
		printf("Error: TexGenome_BuildIndex, inconsistent document number!\n");
		exit(EXIT_FAILURE);
	}

	/* ----- STEP 4 -------- */
	/* match sample to doc index */
	vSampleDocMap = NULL;
	vSampleDocMap = (int *)calloc(nSampleNum, sizeof(int));
	if(vSampleDocMap == NULL)
	{
		printf("Error: TexGenome_BuildIndex, cannot allocate memory for sample-document map!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nSampleNum; ni++)
	{
		vSampleDocMap[ni] = TexGenome_BuildIndex_FindMatch(vGSMIndex[ni]->m_pString, vDocIndex, nDocNum);
		if( (vSampleDocMap[ni] < 0) || (vSampleDocMap[ni] >= nDocNum) )
			printf("Warning: %s has no matching document in the database!\n", vGSMIndex[ni]->m_pString);
	}


	/* release memory */
	for(ni=0; ni<nDocNum; ni++)
	{
		DeleteString(vDocIndex[ni]);
		vDocIndex[ni] = NULL;
	}
	free(vDocIndex);

	for(ni=0; ni<nSampleNum; ni++)
	{
		DeleteString(vGSMIndex[ni]);
		vGSMIndex[ni] = NULL;
	}
	free(vGSMIndex);

	
	/* ----- STEP 5 -------- */
	/* create word-sample score */

	/* allocate memory */
	nLineLength = nDocNum*10;
	vString = NULL;
	vString = (char *)calloc(nLineLength+1, sizeof(char));
	if(vString == NULL)
	{
		printf("Error: TexGenome_BuildIndex, cannot allocate memory for loading word mapping matrix!\n");
		exit(EXIT_FAILURE);
	}

	vDocScore = NULL;
	vDocScore = (float *)calloc(nDocNum, sizeof(float));
	if(vDocScore == NULL)
	{
		printf("Error: TexGenome_BuildIndex, cannot allocate memory for document score!\n");
		exit(EXIT_FAILURE);
	}

	vSampleScore = NULL;
	vSampleScore = (float *)calloc(nSampleNum, sizeof(float));
	if(vSampleScore == NULL)
	{
		printf("Error: TexGenome_BuildIndex, cannot allocate memory for sample score!\n");
		exit(EXIT_FAILURE);
	}

	/* read */
	fpIn = NULL;
	fpIn = fopen(strWordMapFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: TexGenome_BuildIndex, cannot open input word map!\n");
		exit(EXIT_FAILURE);
	}

	nj = 0;
	while(fgets(vString, nLineLength, fpIn) != NULL)
	{
		StrTrimLeft(vString);
		StrTrimRight(vString);
		if(vString[0] == '\0')
			continue;

		if(nj%100 == 0)
			printf(" %d ...\n", nj);

		/* load data */
		chp2 = strpbrk(vString, separators);
		if(chp2 != NULL)
			*chp2 = '\0';
		strcpy(strWord, vString);
		chp1 = chp2+1;
		StrTrimLeft(chp1);

		ni = 0;
		chp2 = strpbrk(chp1, separators);
		while(chp2 != NULL)
		{
			*chp2 = '\0';
			vDocScore[ni] = (float)atof(chp1);
			ni++;

			chp1 = chp2+1;
			StrTrimLeft(chp1);
			chp2 = strpbrk(chp1, separators);
		}
		vDocScore[ni] = (float)atof(chp1);
		ni++;

		if(ni != nDocNum)
		{
			printf("Warning: TexGenome_BuildIndex, inconsistent document number! (%s: %d <> %d)\n", strWord, ni, nDocNum);
			nj++;
			continue;
		}

		/* compute score */
		for(nk=0; nk<nSampleNum; nk++)
		{
			nl = vSampleDocMap[nk];
			if( (nl >= 0) && (nl < nDocNum) )
				vSampleScore[nk] = vDocScore[nl];
			else
				vSampleScore[nk] = 0.0;
		}

		/* write data to file */
		/* printf("%s\t%d\n", strWord, ni); */
		if( ((strWord[0] >= '0') && (strWord[0] <= '9')) || ((strWord[0] >= 'a') && (strWord[0] <= 'z')) )
		{
			if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
				sprintf(strFileName, "%s%c\\%s", strExportFolder, strWord[0], strWord); 
			else
				sprintf(strFileName, "%s%c/%s", strExportFolder, strWord[0], strWord);
		}
		else
		{
			if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
				sprintf(strFileName, "%sother\\%s", strExportFolder, strWord); 
			else
				sprintf(strFileName, "%sother/%s", strExportFolder, strWord);
		}

		TileMapv2_SaveToBinaryFile(vSampleScore, sizeof(float), nSampleNum, strFileName);


		/* prepare next step */
		nj++;
	}

	fclose(fpIn);

	/* ----- STEP 6 -------- */
	/* release memory */
	free(vString);
	free(vDocScore);
	free(vSampleScore);
	free(vSampleDocMap);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TexGenome_BuildIndex_FindMatch()                                       */
/*  Find matching sample in the document index and return the index.       */
/* ----------------------------------------------------------------------- */ 
int TexGenome_BuildIndex_FindMatch(char *strName, struct tagString **vDocIndex, 
				int nDocNum)
{
	/* define */
	int ni,nj,nk;
	int nId = -1;
	int nCmp;
	int nFound = 0;

	ni = 0; 
	nj = nDocNum-1;
	while( (nj-ni)>1 )
	{
		nk = (ni+nj)/2;
		nCmp = strcmp(strName, vDocIndex[nk]->m_pString);

		if(nCmp == 0)
		{
			nId = nk;
			nFound = 1;
			break;
		}
		else if(nCmp < 0)
		{
			nj = nk;
		}
		else
		{
			ni = nk;
		}
	}

	if(nFound == 0)
	{
		if(strcmp(strName, vDocIndex[ni]->m_pString) == 0)
		{
			nId = ni;
			nFound = 1;
		}
		else if(strcmp(strName, vDocIndex[nj]->m_pString) == 0)
		{
			nId = nj;
			nFound = 1;
		}
		else
		{
			nId = -1;
		}
	}

	/* return */
	return nId;
}

/* ----------------------------------------------------------------------- */ 
/*  TexGenome_BuildKeywordProfile()                                        */
/*  ScoreType 0: add; 1: multiply.                                         */
/* ----------------------------------------------------------------------- */ 
int TexGenome_BuildKeywordProfile(char strSampleIndexFile[], char strDictionaryFolder[], 
					char strKeywordFile[], char strOutputFile[], int nScoreType)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	char strKeyword[MED_LINE_LENGTH];
	char *chp1,*chp2;
	char separators[]="\t ";
	int nSampleNum;
	double *vScore;
	int ni,nj;

	/* init */
	AdjustDirectoryPath(strDictionaryFolder);

	fpIn = NULL;
	fpIn = fopen(strSampleIndexFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: TexGenome_BuildKeywordProfile, cannot open sample index file!\n");
		exit(EXIT_FAILURE);
	}

	nSampleNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nSampleNum++;
	}
	fclose(fpIn);

	if(nSampleNum <= 0)
	{
		printf("Error: TexGenome_BuildKeywordProfile, no sample available!");
		return PROC_SUCCESS;
	}

	vScore = NULL;
	vScore = (double *)calloc(nSampleNum, sizeof(double));
	if(vScore == NULL)
	{
		printf("Error: TexGenome_BuildKeywordProfile, cannot create score matrix!");
		return PROC_SUCCESS;
	}

	
	/* open files */
	fpIn = NULL;
	fpIn = fopen(strKeywordFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: TexGenome_BuildIndex, cannot open keyword file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: TexGenome_BuildIndex, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* add score */
		ni = 0;
		chp1 = strLine;
		chp2 = strpbrk(chp1, separators);
		while(chp2 != NULL)
		{
			*chp2 = '\0';
			strcpy(strKeyword, chp1);
			TexGenome_BuildKeywordProfile_ComputeScore(strKeyword,
				strDictionaryFolder, nScoreType, vScore, nSampleNum, ni);

			fprintf(fpOut, "[%s]", strKeyword);
			ni++;
			chp1 = chp2+1;
			StrTrimLeft(chp1);
			chp2 = strpbrk(chp1, separators);
		}

		strcpy(strKeyword, chp1);
		TexGenome_BuildKeywordProfile_ComputeScore(strKeyword,
			strDictionaryFolder, nScoreType, vScore, nSampleNum, ni);
		fprintf(fpOut, "[%s]", strKeyword);
		ni++;

		for(nj=0; nj<nSampleNum; nj++)
		{
			fprintf(fpOut, "\t%e", vScore[nj]);
		}
		fprintf(fpOut, "\n");
	}

	/* close files */
	fclose(fpIn);
	fclose(fpOut);

	/* release memory */
	free(vScore);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TexGenome_BuildKeywordProfile_ComputeScore()                           */
/* ----------------------------------------------------------------------- */ 
int TexGenome_BuildKeywordProfile_ComputeScore(char strKeyword[],
		char strDictionaryFolder[], int nScoreType, double *vScore, 
		int nSampleNum, int nKeyId)
{
	/* define */
	float *vS;
	int nj;
	char strFileName[MED_LINE_LENGTH];
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	FILE *fpIn;

	/* init */
	if(strKeyword[0] == '\0')
	{
		if(nKeyId == 0)
		{
			for(nj=0; nj<nSampleNum; nj++)
				vScore[nj] = 0.0;
		}

		return PROC_SUCCESS;
	}

	vS = NULL;
	vS = (float *)calloc(nSampleNum, sizeof(float));
	if(vS == NULL)
	{
		printf("Error: TexGenome_BuildKeywordProfile_ComputeScore, cannot create vector for word profile!");
		exit(EXIT_FAILURE);
	}

	/* load score */
	if( ((strKeyword[0] >= '0') && (strKeyword[0] <= '9')) || ((strKeyword[0] >= 'a') && (strKeyword[0] <= 'z')) )
	{
		if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
			sprintf(strFileName, "%s%c\\%s", strDictionaryFolder, strKeyword[0], strKeyword); 
		else
			sprintf(strFileName, "%s%c/%s", strDictionaryFolder, strKeyword[0], strKeyword);
	}
	else
	{
		if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
			sprintf(strFileName, "%sother\\%s", strDictionaryFolder, strKeyword); 
		else
			sprintf(strFileName, "%sother/%s", strDictionaryFolder, strKeyword);
	}

	fpIn = NULL;
	fpIn = fopen(strFileName, "rb");
	if(fpIn == NULL)
	{
		printf("Warning: TexGenome_BuildKeywordProfile_ComputeScore, no word [%] in the dictionary!\n");
	}
	else
	{
		if( little_endian_fread(vS, sizeof(float), nSampleNum, fpIn, little_endian_machine) != nSampleNum)
		{
			printf("Error: TexGenome_BuildKeywordProfile_ComputeScore, cannot load data from binary files correctly!\n");
			exit(EXIT_FAILURE);
		}
	}
	fclose(fpIn);

	/* calculate score */
	if(nKeyId == 0)
	{
		for(nj=0; nj<nSampleNum; nj++)
			vScore[nj] = vS[nj];
	}
	else
	{
		for(nj=0; nj<nSampleNum; nj++)
		{
			if(nScoreType == 0)
				vScore[nj] += vS[nj];
			else
				vScore[nj] *= vS[nj];
		}
	}

	/* release memory */
	free(vS);
	
	/* return */
	return PROC_SUCCESS;
}