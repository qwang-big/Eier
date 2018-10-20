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
/*  TexGenome_BuildIndex()                                                 */
/* ----------------------------------------------------------------------- */ 
int TexGenome_BuildIndex(char strSampleIndexFile[], char strDocIndexFile[], 
						 char strWordMapFile[], char strExportFolder[]);

/* ----------------------------------------------------------------------- */ 
/*  TexGenome_BuildIndex_FindMatch()                                       */
/*  Find matching sample in the document index and return the index.       */
/* ----------------------------------------------------------------------- */ 
int TexGenome_BuildIndex_FindMatch(char *strName, struct tagString **vDocIndex, 
				int nDocNum);

/* ----------------------------------------------------------------------- */ 
/*  TexGenome_BuildKeywordProfile()                                        */
/*  ScoreType 0: add; 1: multiply.                                         */
/* ----------------------------------------------------------------------- */ 
int TexGenome_BuildKeywordProfile(char strSampleIndexFile[], char strDictionaryFolder[], 
					char strKeywordFile[], char strOutputFile[], int nScoreType);

/* ----------------------------------------------------------------------- */ 
/*  TexGenome_BuildKeywordProfile_ComputeScore()                           */
/* ----------------------------------------------------------------------- */ 
int TexGenome_BuildKeywordProfile_ComputeScore(char strKeyword[],
		char strDictionaryFolder[], int nScoreType, double *vScore, 
		int nSampleNum, int nKeyId);

#ifdef __cplusplus
}
#endif
