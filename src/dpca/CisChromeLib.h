#ifdef __cplusplus
extern "C" {
#endif
/* ----------------------------------------------------------------------- */
/*                                 Macro                                   */
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */ 
/*                                 Struct                                  */
/* ----------------------------------------------------------------------- */ 


/* ----------------------------------------------------------------------- */
/*                              Declaration                                */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_Main()                                                       */
/*  CisChrome: analyzing chromatin signals at DNA motif sites.             */
/* ----------------------------------------------------------------------- */ 
int CisChrome_Main(char strParamFile[]);

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
		int *pGridNum, double *pFDRCut, int *pTSS, int *pUp, int *pDown);

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_ReadTFInfo()                                                 */
/*  Read TF information.                                                   */
/* ----------------------------------------------------------------------- */
int CisChrome_ReadTFInfo(char strMotifLibrary[], int *pMotifNum, 
	struct tagString ***pvTFName, struct tagRefGene ***vTFRefGene,
	struct tagString ***pvMotifFilePath, int **pvTFEntrezID, char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_GenomeProfile()                                              */
/*  Create Genome-wide enrichment profiles.                                */
/* ----------------------------------------------------------------------- */
int CisChrome_GenomeProfile(int nGroupNum, int nDatasetNum, int nSampleNum, 
		int nPaired, int nTwosided, int *vGroupId, int *vDatasetId, int *vExtLen, 
		int **vRepNum, struct tagString **vSampleName, struct tagString **vSampleFilePath,
		int nBinSize, int nWinSize, char strWorkFolder[], char strProjectTitle[],
		int *pChrNum, struct tagString ***vChrName, struct INTMATRIX **pChrLen,
		int nGridNum, double *pMinT, double *pMaxT);

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_GenomeProfile_Initialize()                                   */
/*  Initialize genomewide enrichment profiles.                             */
/* ----------------------------------------------------------------------- */
int CisChrome_GenomeProfile_Initialize(int nSampleNum, 
		struct tagString **vSampleName, struct tagString **vSampleFilePath, 
		int *vExtLen, int nBinSize, struct DOUBLEMATRIX **ppFileReadCount, 
		struct DOUBLEMATRIX **ppFileNormFactor, 
		int *pBaseReadNum, int *pBaseLineFileId, 
		int *pChrNum, struct tagString ***vvChrName, struct INTMATRIX **ppChrLen);

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_GenomeProfile_CountBinReads()                                */
/*  Count reads for genomic bins                                           */
/* ----------------------------------------------------------------------- */ 
int CisChrome_GenomeProfile_CountBinReads(int nSampleNum, 
		struct tagString **vSampleName, struct tagString **vSampleFilePath, 
		int *vExtLen, int nBinSize,	struct DOUBLEMATRIX *pFileNormFactor,
		char strWorkFolder[], char strProjectTitle[],  
		int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen,
		int nExportBAR);

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_GenomeProfile_CountBinReads_SingleFile()                     */
/*  Count reads for genomic bins for a single sample file                  */
/* ----------------------------------------------------------------------- */ 
int CisChrome_GenomeProfile_CountBinReads_SingleFile(char strInputFile[], char strOutputPath[],
				char strOutputFile[], int nChrNum, struct tagString **vChrName, 
				struct INTMATRIX *pChrLen, int nBinSize, 		
				int nExtLen, double dNormFactor, int nExportBAR);

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_ComputeBinStats()                                            */
/*  Compute summary statistics for genomic bins.                           */
/* ----------------------------------------------------------------------- */ 
int CisChrome_ComputeBinStats(int nGroupNum, int nDatasetNum, int nSampleNum, 
		int nPaired, int nTwosided, int *vGroupId, int *vDatasetId, int **vRepNum,
		struct tagString **vSampleName, int nBinSize, int nWinSize,
		char strWorkFolder[], char strProjectTitle[], 
		int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen,
		int nGridNum, double *pMinT, double *pMaxT);

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_ComputeBinStats_Chr_Initial()                                */
/*  Compute initial summary statistics for genomic bins for a single       */
/*  chromosome.                                                            */
/* ----------------------------------------------------------------------- */ 
int CisChrome_ComputeBinStats_Chr_Initial(int nGroupNum, int nDatasetId,
			int nSampleNum, int nPaired, int *vGroupId, int *vDatasetId, 
			int **vRepNum, struct tagString **vSampleName, int nBinSize, 
			int nWinSize, int nBinNum, char strWorkFolder[], char strProjectTitle[],
			char strChr[], double *pS2M, double *pS2S, int *pEffectBinNum);

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_ComputeBinStats_Chr_Tstat()                                  */
/*  Compute t-statistics.                                                  */
/* ----------------------------------------------------------------------- */ 
int CisChrome_ComputeBinStats_Chr_Tstat(int nDatasetId, char strWorkFolder[], 
			char strProjectTitle[], char strChr[], int nBinNum, double dB, 
			double dVM, double dVScale, double *pMinT, double *pMaxT);

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_ComputeBinStats_Chr_Thist()                                  */
/*  Get distribution of t-statistics.                                      */
/* ----------------------------------------------------------------------- */ 
int CisChrome_ComputeBinStats_Chr_Thist(int nDatasetId, char strWorkFolder[], 
		char strProjectTitle[],  char strChr[], int nBinNum, int nGridNum,
		double dMinT, double dMaxT, double dStep, float *vThist);

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
		struct tagRefGene **vTFRefGene, int nUp, int nDown, int nTSS);

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_FindActiveTFBS_LoadTFBS()                                    */
/*  Load motif sites.                                                      */
/* ----------------------------------------------------------------------- */ 
int CisChrome_FindActiveTFBS_LoadTFBS(char strMotifFileName[], int *pTFBSNum, 
		struct tagString ***pvTFBSChr, int **pvTFBSPos1, int **pvTFBSPos2,
		unsigned char **pvTFBSStrand);

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_FindActiveTFBS_Call()                                        */
/*  Call active motif sites.                                               */
/* ----------------------------------------------------------------------- */ 
int CisChrome_FindActiveTFBS_Call(char strWorkFolder[], char strProjectTitle[], 
		int nDatasetId, int nTFId, char strTFName[], int nTFBSNum, 
		struct tagString **vTFBSChr, int *vTFBSPos1, int *vTFBSPos2, 
		unsigned char *vTFBSStrand, int nChrNum, struct tagString **vChrName, 
		struct INTMATRIX *pChrLen, int nBinSize, int nWinSize, int nTwoSided,
		float *vMotifHist, float *vMotifCDF, double dMinT, double dMaxT, 
		double dStepSize, int nGridNum, double dFDRCut,
		char strTFChr[], int nTFPos1, int nTFPos2, double *pTFExpress);

/* ----------------------------------------------------------------------- */ 
/*  CisChrome_Export()                                                     */
/*  Export results.                                                        */
/* ----------------------------------------------------------------------- */ 
int CisChrome_Export(char strWorkFolder[], char strProjectTitle[],
		int nMotifNum, int nDatasetNum, int nGridNum, struct tagString **vTFName, 
		int *vTotalSiteNum, int *vActiveSiteNum, double *vActiveSitePrc,
		float **vMotifStat, double *vMinT, double *vMaxT, 
		double *vTFExpress, int nOutputType);
#ifdef __cplusplus
}
#endif
