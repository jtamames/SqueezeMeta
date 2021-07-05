/*System includes*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <float.h>

/*GSL includes*/
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <pthread.h>
#include <gsl/gsl_sf_exp.h>

/*User includes*/
#include "c_vbgmm_fit.h"

void readInputData(const char *szFile, t_Data *ptData);

void c_vbgmm_fit (double* adX, int nN, int nD, int nK, int* anAssign, int debug, int bAssign);

void readAssigns(const char *szFile, int *anAssign, int nN);

int main()
{
    t_Data tData;
    int *anAssign = NULL;
    int i = 0, j = 0, nN = -1, nD = -1, nK = 0;
    double* adX = NULL;
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);    
    
    
    readInputData("PCA_transformed_data_gt2000.csv", &tData);
    
    nN = tData.nN;
    nD = tData.nD;
    
    adX = (double *) malloc(nN*nD*sizeof(double));
    anAssign = (int *) malloc(nN*sizeof(int));
    
    readAssigns("clustering_gt2000.csv", anAssign, nN);
    
    for(i = 0; i < nN; i++){
        if(anAssign[i] > nK){
            nK = anAssign[i];
        }
        for(j = 0; j < nD; j++){
            adX[i*nD + j] = tData.aadX[i][j];    
        }
    }
    
    nK = nK + 1;
    fprintf(stderr,"Run c_vbgmm_fit with %d clusters\n",nK);
    fflush(stderr);
    
    c_vbgmm_fit (adX, nN, nD, nK, anAssign, FALSE, TRUE);
    for(i = 0; i < nN; i++){
        printf("%d,%d\n",i,anAssign[i]);
    }
    free(adX);
    free(anAssign);
    return 0;
}

void readInputData(const char *szFile, t_Data *ptData)
{
  double  **aadX = NULL;
  int  i = 0, j = 0, nD = 0, nN = 0;
  char *szLine = (char *) malloc(sizeof(char)*MAX_LINE_LENGTH);
  FILE* ifp = NULL;

  if(!szLine)
    goto memoryError;

  ifp = fopen(szFile, "r");

  if(ifp){
    char* szTok   = NULL;
    char* pcError = NULL;

    if(fgets(szLine, MAX_LINE_LENGTH, ifp) == NULL)
      goto formatError;

    szTok = strtok(szLine, DELIM);
    /*count dimensions*/
    while(strtok(NULL, DELIM) != NULL){
      
      nD++;
    }
    /*count data points*/
    while(fgets(szLine, MAX_LINE_LENGTH, ifp) != NULL){
    	nN++;
    }
    fclose(ifp);

    /*reopen input file*/
    ifp = fopen(szFile, "r");	

    if(fgets(szLine, MAX_LINE_LENGTH, ifp) == NULL)
      goto formatError;


    /*allocate memory for dimension names*/
    //ptData->aszDimNames = (char **) malloc(nD*sizeof(char*));
    //if(!ptData->aszDimNames)
      //goto memoryError;

    szTok = strtok(szLine, DELIM);
    /*read in dim names*/
    for(i = 0; i < nD; i++){
      szTok = strtok(NULL, DELIM);
     // ptData->aszDimNames[i] = strdup(szTok);
    }
	
    /*allocate memory for data matrix*/
    aadX = (double **) malloc(nN*sizeof(double*));
    if(!aadX)
      goto memoryError;
    for(i = 0; i < nN; i++){
      aadX[i] = (double *) malloc(nD*sizeof(double));
      if(!aadX[i])
	goto memoryError;
    }

    /*read in input data*/
    //ptData->aszSampleNames = (char **) malloc(nN*sizeof(char*));
    //if(!ptData->aszSampleNames)
      //goto memoryError;

    for(i = 0; i < nN; i++){
    
      if(fgets(szLine, MAX_LINE_LENGTH, ifp) == NULL)
	goto formatError;

      szTok = strtok(szLine, DELIM);
     // ptData->aszSampleNames[i] = strdup(szTok);
      for(j = 0; j < nD; j++){
	szTok = strtok(NULL, DELIM);

	aadX[i][j] = strtod(szTok,&pcError);

	if(*pcError != '\0'){
	  goto formatError;
	}
      }
    }
  }
  else{
    fprintf(stderr, "Failed to open abundance data file %s aborting\n", szFile);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }

  free(szLine);
  ptData->nD = nD;
  ptData->nN = nN;
  ptData->aadX = aadX;
  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in readInputData\n");
  fflush(stderr);
  exit(EXIT_FAILURE);

 formatError:
  fprintf(stderr, "Incorrectly formatted abundance data file\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

void readAssigns(const char *szFile, int *anAssign, int nN)
{
  int  i = 0;
  char *szLine = (char *) malloc(sizeof(char)*MAX_LINE_LENGTH);
  FILE* ifp = NULL;

  if(!szLine)
    goto memoryError;

  ifp = fopen(szFile, "r");

  if(ifp){
    char* szTok   = NULL;
    char* pcError = NULL;

    if(fgets(szLine, MAX_LINE_LENGTH, ifp) == NULL)
        goto formatError;

    for(i = 0; i < nN; i++){ 
        if(fgets(szLine, MAX_LINE_LENGTH, ifp) == NULL)
            goto formatError;

        szTok = strtok(szLine, DELIM);
        szTok = strtok(NULL, DELIM);
    	anAssign[i] = strtod(szTok,&pcError);

	    if(*pcError != '\0'){
	        goto formatError;
	    }
    }
  }
  else{
    fprintf(stderr, "Failed to open abundance data file %s aborting\n", szFile);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }

  free(szLine);

  return;

 memoryError:
  fprintf(stderr, "Failed allocating memory in readInputData\n");
  fflush(stderr);
  exit(EXIT_FAILURE);

 formatError:
  fprintf(stderr, "Incorrectly formatted abundance data file\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}
