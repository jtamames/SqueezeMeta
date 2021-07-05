/* C functions for running vbgmm from Cython*/

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
#include <omp.h>

/*User includes*/
#include "c_vbgmm_fit.h"

void c_vbgmm_fit (double* adX, int nN, int nD, int nK, int seed, int* anAssign, int nThreads, int nIter)
{
    int debug = 0;
    int bAssign = 0;

    if (nIter < 1){
        nIter = DEF_MAX_ITER;
    }

    driverMP(adX, nN, nD, anAssign, nK, seed, nIter, DEF_EPSILON, debug, bAssign, nThreads);

    return;
}

int driverMP(double *adX, int nN, int nD, int *anAssign, int nKStart, unsigned long lSeed, 
                                        int nMaxIter, double dEpsilon, int debug, int bAssign, int nThreads)
{
    t_Params           tParams;
    t_Data             tData;
    gsl_rng            *ptGSLRNG     = NULL;
    const gsl_rng_type *ptGSLRNGType = NULL;
    t_VBParams tVBParams;
    t_Cluster  *ptCluster = NULL;
    int i = 0, nK = nKStart, nNthreads = 0;
    int nT = 0;
    char *szCOutFile = NULL;

    /*initialise GSL RNG*/
    gsl_rng_env_setup();

    gsl_set_error_handler_off();

    ptGSLRNGType = gsl_rng_default;
    ptGSLRNG     = gsl_rng_alloc(ptGSLRNGType);

    /*set OMP thread number*/
    nNthreads = nThreads;
    nT = nN / 32  + 1;
    printf("%d %d %d\n",nN,nT,nNthreads);
    if (nT < nNthreads){
        nNthreads = nT;
    }
    
    omp_set_num_threads(nNthreads);
    fprintf(stderr,"Setting %d OMP threads\n",nNthreads);

    /*set clusters params*/
    tParams.nKStart  = nKStart;
    tParams.nMaxIter = nMaxIter;
    tParams.dEpsilon = dEpsilon;
    tParams.lSeed    = lSeed;
    
    fprintf(stderr,"Generate input data\n");
    fflush(stderr);
    generateInputData(adX, nN, nD, &tData);

    setVBParams(&tVBParams, &tData);

    if(debug>0){
        szCOutFile = (char *) malloc(sizeof(char)*MAX_FILE_NAME_LENGTH);
        
	    sprintf(szCOutFile,"%sr%d.csv",DEF_FILE_STUB,0);
	}

    if(bAssign > 0){
        nK = 0;
        
        for(i = 0; i < nN; i++){
            if(anAssign[i] > nK){
                nK = anAssign[i];
            }
        }
    }
    ptCluster = malloc(sizeof(t_Cluster));
    allocateCluster(ptCluster,nN,nK,nD,&tData,lSeed,nMaxIter,dEpsilon,szCOutFile);

    ptCluster->bAssign = bAssign;
    if(bAssign > 0){        
        for(i = 0; i < nN; i++){
            ptCluster->anMaxZ[i] = anAssign[i];
        }
    }

    ptCluster->ptVBParams = &tVBParams;
        
    ptCluster->nThread = 0;

    fitEM_MP((void *) ptCluster);

    compressCluster(ptCluster);

    calcCovarMatrices(ptCluster,&tData);

    for(i = 0; i < nN; i++){
      anAssign[i] = ptCluster->anMaxZ[i];
    }

    if(debug>0){
        free(szCOutFile);
	}

    /*free up memory in data object*/
    destroyData(&tData);

    /*free up best BIC clusters*/

    destroyCluster(ptCluster);
    free(ptCluster);

    gsl_rng_free(ptGSLRNG);
    gsl_matrix_free(tVBParams.ptInvW0);

    return EXIT_SUCCESS;
    
    memoryError:
        fprintf(stderr, "Failed allocating memory in driver\n");
        fflush(stderr);
        exit(EXIT_FAILURE);
}

void generateInputData(double *adX, int nN, int nD, t_Data *ptData)
{
    double  **aadX = NULL;
    int  i = 0, j = 0;

    /*allocate memory for data matrix*/
    aadX = (double **) malloc(nN*sizeof(double*));
    if(!aadX)
        goto memoryError;

    for(i = 0; i < nN; i++){
        aadX[i] = (double *) malloc(nD*sizeof(double));
        if(!aadX[i])
	        goto memoryError;
    }

    for(i = 0; i < nN; i++){
        for(j = 0; j < nD; j++){
            aadX[i][j] = adX[i*nD + j];
        }
    }
    ptData->nD = nD;
    ptData->nN = nN;
    ptData->aadX = aadX;
    return;

 memoryError:
    fprintf(stderr, "Failed allocating memory in readInputData\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
}

void destroyData(t_Data *ptData)
{
    int nN = ptData->nN;
    int i = 0;

    for(i = 0; i < nN; i++){
        free(ptData->aadX[i]);
    }
    free(ptData->aadX);

    return;
}

void calcSampleVar(t_Data *ptData,double *adVar, double *adMu)
{
    double **aadX = ptData->aadX;
    int i = 0, n = 0;
    int nD = ptData->nD, nN = ptData->nN;
    /*sample means*/
    double dN = (double) nN;

    for(i = 0; i < nD; i++){
        adMu[i] = 0.0;
        adVar[i] = 0.0;
    }

    for(i = 0; i < nD; i++){
        for(n = 0; n < nN; n++){
            adMu[i] += aadX[n][i];
            adVar[i] += aadX[n][i]*aadX[n][i];
        }

        adMu[i] /= dN;

        adVar[i] = (adVar[i] - dN*adMu[i]*adMu[i])/(dN - 1.0);
    }

    return;
}

void setVBParams(t_VBParams *ptVBParams, t_Data *ptData)
{
    int i = 0, nD = ptData->nD;
    double adVar[nD], adMu[nD];

    ptVBParams->dBeta0 = DEF_BETA0;
    ptVBParams->dNu0 = (double) nD;
    ptVBParams->ptInvW0 = gsl_matrix_alloc(nD,nD);

    calcSampleVar(ptData,adVar, adMu);
    gsl_matrix_set_zero (ptVBParams->ptInvW0);

    for(i = 0; i < nD; i++){
        double dRD = adVar[i]*((double) nD);

        gsl_matrix_set(ptVBParams->ptInvW0,i,i,dRD);
    }

    ptVBParams->dLogWishartB = dLogWishartB(ptVBParams->ptInvW0, nD, ptVBParams->dNu0, TRUE);
}

void allocateCluster(t_Cluster *ptCluster, int nN, int nK, int nD, t_Data *ptData, long lSeed, int nMaxIter, double dEpsilon, char *szCOutFile)
{
    int i = 0, j = 0, k = 0;

    ptCluster->szCOutFile = szCOutFile;
    ptCluster->ptVBParams = NULL;
    ptCluster->lSeed = lSeed;
    ptCluster->nMaxIter = nMaxIter;
    ptCluster->dEpsilon = dEpsilon;
    ptCluster->ptData = ptData;

    ptCluster->nN = nN;
    ptCluster->nK = nK;
    ptCluster->nKSize = nK;
    ptCluster->nD = nD;

    ptCluster->dVBL = 0.0;

    ptCluster->anMaxZ = (int *) malloc(nN*sizeof(int)); /*destroyed*/
    if(!ptCluster->anMaxZ)
        goto memoryError;

    ptCluster->anW = (int *) malloc(nK*sizeof(int)); /*destroyed*/
    if(!ptCluster->anW)
        goto memoryError;

    for(i = 0; i < nN; i++){
        ptCluster->anMaxZ[i] = NOT_SET;
    }

    for(i = 0; i < nK; i++){
        ptCluster->anW[i] = 0;
    }

    ptCluster->aadZ = (double **) malloc(nN*sizeof(double *)); /*destroyed*/
    if(!ptCluster->aadZ)
        goto memoryError;

    for(i = 0; i < nN; i++){
        ptCluster->aadZ[i] = (double *) malloc(nK*sizeof(double)); /*destroyed*/
        if(!ptCluster->aadZ[i])
            goto memoryError;

        for(j = 0; j < nK; j++){
            ptCluster->aadZ[i][j] = 0.0;
        }
    }

    ptCluster->adLDet = (double *) malloc(nK*sizeof(double)); /*all*/
    ptCluster->adPi   = (double *) malloc(nK*sizeof(double));
    ptCluster->adBeta = (double *) malloc(nK*sizeof(double));
    ptCluster->adNu   = (double *) malloc(nK*sizeof(double)); /*destroyed*/

    if(!ptCluster->adLDet || !ptCluster->adPi)
        goto memoryError;

    if(!ptCluster->adBeta || !ptCluster->adNu)
        goto memoryError;

    for(k = 0; k < nK; k++){
        ptCluster->adLDet[k] = 0.0;
        ptCluster->adPi[k]   = 0.0;
        ptCluster->adBeta[k] = 0.0;
        ptCluster->adNu[k] = 0.0;
    }

    ptCluster->aadMu = (double **) malloc(nK*sizeof(double *));
    if(!ptCluster->aadMu)
        goto memoryError;

    ptCluster->aadM = (double **) malloc(nK*sizeof(double *));
    if(!ptCluster->aadM)
        goto memoryError;

    for(i = 0; i < nK; i++){
        ptCluster->aadM[i] = (double*) malloc (nD*sizeof(double));
        if(!ptCluster->aadM[i])
            goto memoryError;

        ptCluster->aadMu[i] = (double*) malloc (nD*sizeof(double));
        if(!ptCluster->aadMu[i])
            goto memoryError;
    }

    ptCluster->aptSigma = (gsl_matrix **) malloc(nK*sizeof(gsl_matrix *));
    if(!ptCluster->aptSigma)
        goto memoryError;

    for(i = 0; i < nK ; i++){
        ptCluster->aptSigma[i] = (gsl_matrix*) gsl_matrix_alloc (nD, nD);
    }

    ptCluster->aptCovar = (gsl_matrix **) malloc(nK*sizeof(gsl_matrix *));
    if(!ptCluster->aptCovar)
        goto memoryError;

    for(i = 0; i < nK ; i++){
        ptCluster->aptCovar[i] = (gsl_matrix*) gsl_matrix_alloc (nD, nD);
    }

    return;

    memoryError:
    fprintf(stderr, "Failed allocating memory in allocateCluster\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
}

void destroyCluster(t_Cluster* ptCluster)
{
  int i = 0, nN = ptCluster->nN, nKSize = ptCluster->nKSize;

  if(ptCluster->szCOutFile != NULL){
	  free(ptCluster->szCOutFile);
  }

  free(ptCluster->anMaxZ);

  free(ptCluster->anW);

  for(i = 0; i < nN; i++){
    free(ptCluster->aadZ[i]);
  }
  free(ptCluster->aadZ);

  free(ptCluster->adLDet);
  free(ptCluster->adPi);
  free(ptCluster->adBeta);
  free(ptCluster->adNu);

  for(i = 0; i < nKSize; i++){
    free(ptCluster->aadMu[i]);
    free(ptCluster->aadM[i]);
  }

  free(ptCluster->aadMu);
  free(ptCluster->aadM);

  for(i = 0; i < nKSize ; i++){
    gsl_matrix_free(ptCluster->aptSigma[i]);
    gsl_matrix_free(ptCluster->aptCovar[i]);
  }
  free(ptCluster->aptSigma);
  free(ptCluster->aptCovar);
  return;
}

void* fitEM_MP(void *pvCluster)
{
    t_Cluster          *ptCluster = (t_Cluster *) pvCluster;
    gsl_rng            *ptGSLRNG     = NULL;
    const gsl_rng_type *ptGSLRNGType = NULL;
    int i = 0, k = 0;
    /*initialise GSL RNG*/
    ptGSLRNGType = gsl_rng_default;
    ptGSLRNG     = gsl_rng_alloc(ptGSLRNGType);

    gsl_rng_set (ptGSLRNG, ptCluster->lSeed);

    if(ptCluster->bAssign == FALSE){
        initKMeans(ptGSLRNG, ptCluster, ptCluster->ptData);
    }
    else{
        for(i = 0; i < ptCluster->nN; i++){
            for(k = 0; k < ptCluster->nK; k++){
                ptCluster->aadZ[i][k] = 0.0;
            }
            ptCluster->aadZ[i][ptCluster->anMaxZ[i]] = 1.0;
        }
        performMStepMP(ptCluster, ptCluster->ptData);
    }
    
    gmmTrainVB_MP(ptCluster, ptCluster->ptData);

    gsl_rng_free(ptGSLRNG);

    return NULL;
}

void compressCluster(t_Cluster *ptCluster)
{
    int i = 0, k = 0, nNewK = 0, nN = ptCluster->nN;
    double **aadNewZ = NULL, dN = (double) nN;

    for(i = 0; i < ptCluster->nK; i++){
        if(ptCluster->adPi[i] > 0.0){
            nNewK++;
        }
    }

    aadNewZ = (double **) malloc(nN*sizeof(double *));
    if(!aadNewZ)
        goto memoryError;

    for(i = 0; i < nN; i++){
        aadNewZ[i] = (double *) malloc(nNewK*sizeof(double));
        if(!aadNewZ[i])
            goto memoryError;
    }

    for(i = 0; i < nN; i++){
        int nC = 0;
        for(k = 0; k < ptCluster->nK; k++){
            if(ptCluster->adPi[k] > 0.0){
	            aadNewZ[i][nC] = ptCluster->aadZ[i][k];
	            nC++;
            }
        }
    }

    for(i = 0; i < nN; i++){
        free(ptCluster->aadZ[i]);
    }
    free(ptCluster->aadZ);

    /*reset Z and K*/
    ptCluster->aadZ = aadNewZ;
    ptCluster->nK = nNewK;

    /*recalculate Pi*/
    for(k = 0; k < ptCluster->nK; k++){
        ptCluster->adPi[k] = 0.0;
        for(i = 0; i < nN; i++){
            ptCluster->adPi[k] += ptCluster->aadZ[i][k];
        }
        ptCluster->adPi[k] /= dN;
    }

    /*assign to best clusters*/
    for(i = 0; i < nN; i++){
        double dMaxZ = ptCluster->aadZ[i][0];
        int    nMaxK = 0;
        for(k = 1; k < ptCluster->nK; k++){
            if(ptCluster->aadZ[i][k] > dMaxZ){
	            nMaxK = k;
	            dMaxZ = ptCluster->aadZ[i][k];
            }
        }
        ptCluster->anMaxZ[i] = nMaxK;
    }

    for(k = 0; k < ptCluster->nK; k++){
        ptCluster->anW[k] = 0;
    }

    for(i = 0; i < nN; i++){
        ptCluster->anW[ptCluster->anMaxZ[i]]++;
    }

    return;

    memoryError:
    fprintf(stderr, "Failed allocating memory in main\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
}

double decomposeMatrix(gsl_matrix *ptSigmaMatrix, int nD)
{
  double dDet = 0.0;
  int status;
  int l = 0;

  status = gsl_linalg_cholesky_decomp(ptSigmaMatrix);

  if(status == GSL_EDOM){
    fprintf(stderr,"Failed Cholesky decomposition in decomposeMatrix\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  else{
    for(l = 0; l < nD; l++){
      double dT = gsl_matrix_get(ptSigmaMatrix,l,l);
      dDet += 2.0*log(dT);
    }
    gsl_linalg_cholesky_invert(ptSigmaMatrix);
    return dDet;
  }
}

double mstep(int k, double *adMu, double* adM, int nN, int nD, double** aadZ, double** aadX, double *pdBeta, double *pdNu, double *pdLDet, t_VBParams *ptVBParams, gsl_matrix *ptCovarK, gsl_matrix *ptSigmaMatrix)
{
    double dPi = 0.0, dBeta = 0.0, dLDet = 0.0, dNu = 0.0;
    int i = 0, j = 0, l = 0, m = 0;
    double      dF = 0.0;
    double **aadCovar = NULL;
    double **aadInvWK = NULL;
    
    aadCovar = (double **) malloc(nD*sizeof(double*));
    if(!aadCovar)
        goto memoryError;

    for(i = 0; i < nD; i++){
        aadCovar[i] = (double *) malloc(nD*sizeof(double));
        if(!aadCovar[i])
            goto memoryError;
    }
    
    aadInvWK = (double **) malloc(nD*sizeof(double*));
    if(!aadInvWK)
        goto memoryError;

    for(i = 0; i < nD; i++){
        aadInvWK[i] = (double *) malloc(nD*sizeof(double));
        if(!aadInvWK[i])
            goto memoryError;
    }
    
    /*recompute mixture weights and means*/
    for(j = 0; j < nD; j++){
        adMu[j] = 0.0;
        for(l = 0; l < nD; l++){
            aadCovar[j][l] = 0.0;
            aadInvWK[j][l] = 0.0;
        }
    }

    for(i = 0; i < nN; i++){
        if(aadZ[i][k] > MIN_Z){
            dPi += aadZ[i][k];
            for(j = 0; j < nD; j++){
                adMu[j] += aadZ[i][k]*aadX[i][j];
            }
        }
    }

    /*normalise means*/
    if(dPi > MIN_PI){
        /*Equation 10.60*/
        dBeta = ptVBParams->dBeta0 + dPi;
        
        for(j = 0; j < nD; j++){
            /*Equation 10.61*/
            adM[j] = adMu[j]/dBeta;
            adMu[j] /= dPi;
        }

        dNu = ptVBParams->dNu0 + dPi;

        /*calculate covariance matrices*/
        for(i = 0; i < nN; i++){
            if(aadZ[i][k] > MIN_Z){
                double adDiff[nD];

                for(j = 0; j < nD; j++){
                    adDiff[j] = aadX[i][j] - adMu[j];
                }

                for(l = 0; l < nD; l++){
                    for(m = 0; m <=l ; m++){
                        aadCovar[l][m] += aadZ[i][k]*adDiff[l]*adDiff[m];
                    }
                }
            }
        }

        for(l = 0; l < nD; l++){
            for(m = l + 1; m < nD; m++){
                aadCovar[l][m] = aadCovar[m][l];
            }
        }

        /*save sample covariances for later use*/
        for(l = 0; l < nD; l++){
            for(m = 0; m < nD; m++){
                double dC = aadCovar[l][m] / dPi;
                gsl_matrix_set(ptCovarK,l,m,dC);
            }
        }

        /*Now perform equation 10.62*/
        dF = (ptVBParams->dBeta0*dPi)/dBeta;
        for(l = 0; l < nD; l++){
            for(m = 0; m <= l; m++){
                aadInvWK[l][m] = gsl_matrix_get(ptVBParams->ptInvW0, l,m) + aadCovar[l][m] + dF*adMu[l]*adMu[m];
            }
        }

        for(l = 0; l < nD; l++){
            for(m = 0; m <= l ; m++){
                aadCovar[l][m] /= dPi;
                gsl_matrix_set(ptSigmaMatrix, l, m, aadInvWK[l][m]);
                gsl_matrix_set(ptSigmaMatrix, m, l, aadInvWK[l][m]);
            }
       }


    /*Implement Equation 10.65*/
      dLDet = ((double) nD)*log(2.0);

      for(l = 0; l < nD; l++){
        double dX = 0.5*(dNu - (double) l);
        dLDet += gsl_sf_psi (dX);
      }

      dLDet -= decomposeMatrix(ptSigmaMatrix,nD);
    }
    else{
      /*Equation 10.60*/
      dPi = 0.0;

      dBeta = ptVBParams->dBeta0;

      for(j = 0; j < nD; j++){
        /*Equation 10.61*/
        adM[j] = 0.0;
        adMu[j] = 0.0;
      }

      dNu = ptVBParams->dNu0;

      for(l = 0; l < nD; l++){
        for(m = 0; m <= l; m++){
          aadInvWK[l][m] = gsl_matrix_get(ptVBParams->ptInvW0, l,m);
        }
      }

      for(l = 0; l < nD; l++){
        for(m = 0; m <= l ; m++){
            aadInvWK[l][m] = gsl_matrix_get(ptVBParams->ptInvW0, l,m);
        }
      }

      for(l = 0; l < nD; l++){
        for(m = 0; m <= l ; m++){
            gsl_matrix_set(ptSigmaMatrix, l, m, aadInvWK[l][m]);
            gsl_matrix_set(ptSigmaMatrix, m, l, aadInvWK[l][m]);
        }
      }

      /*Implement Equation 10.65*/
      dLDet = ((double) nD)*log(2.0);

      for(l = 0; l < nD; l++){
        double dX = 0.5*(dNu - (double) l);
        dLDet += gsl_sf_psi (dX);
      }

      dLDet -= decomposeMatrix(ptSigmaMatrix,nD);
    }
    
    /*free up memory*/
    for(i = 0; i < nD; i++){
        free(aadCovar[i]);
        free(aadInvWK[i]);
    }

    free(aadCovar);
    free(aadInvWK);
    
    (*pdBeta) = dBeta;
    (*pdNu)   = dNu;
    (*pdLDet) = dLDet;
    return dPi;

    memoryError:
        fprintf(stderr, "Failed allocating memory in mstep\n");
        fflush(stderr);
        exit(EXIT_FAILURE);
}

void performMStepMP(t_Cluster *ptCluster, t_Data *ptData){
    int k = 0;
    int nN = ptData->nN, nK = ptCluster->nK, nD = ptData->nD;
    double **aadZ = ptCluster->aadZ,**aadX = ptData->aadX;
    double **aadCovar = NULL, **aadInvWK = NULL;
    t_VBParams *ptVBParams = ptCluster->ptVBParams;

    /*perform M step*/
#pragma omp parallel for
    for(int k = 0; k < nK; k++){ /*loop components*/
        double dPi = 0.0, dBeta = 0.0, dNu = 0.0, dLDet = 0.0;

        dPi = mstep(k, ptCluster->aadMu[k],ptCluster->aadM[k], nN, nD, aadZ, aadX, &dBeta, &dNu, &dLDet, ptVBParams, ptCluster->aptCovar[k], ptCluster->aptSigma[k]);
        
        ptCluster->adPi[k] = dPi;
        ptCluster->adBeta[k] = dBeta;
        ptCluster->adNu[k] = dNu;
        ptCluster->adLDet[k] = dLDet;
    }

    /*Normalise pi*/

    if(1){
        double dNP = 0.0;

        for(k = 0; k < nK; k++){
            dNP += ptCluster->adPi[k];
        }

        for(k = 0; k < nK; k++){
            ptCluster->adPi[k] /= dNP;
        }
    }

    return;

    memoryError:
        fprintf(stderr, "Failed allocating memory in performMStep\n");
        fflush(stderr);
        exit(EXIT_FAILURE);
}

void initKMeans(gsl_rng *ptGSLRNG, t_Cluster *ptCluster, t_Data *ptData)
{
    /*very simple initialisation assign each data point to random cluster*/
    int i = 0, k = 0, nN = ptData->nN, nK = ptCluster->nK, nD = ptData->nD;
    double **aadMu = ptCluster->aadMu, **aadX = ptData->aadX;
    int *anMaxZ = ptCluster->anMaxZ, *anW = ptCluster->anW, nChange = nN;
    int nIter = 0, nMaxIter = ptCluster->nMaxIter;
    for(i = 0; i < nN; i++){
        int nIK = gsl_rng_uniform_int (ptGSLRNG, nK);
        ptCluster->anMaxZ[i] = nIK;
        anW[nIK]++;
    }

    updateMeans(ptCluster, ptData);

    while(nChange > 0 && nIter < nMaxIter){
        nChange = 0;
        /*reassign vectors*/
        for(i = 0; i < nN; i++){
        double dMinDist = DBL_MAX;
        int    nMinK = NOT_SET;

        for(k = 0; k < nK; k++){
	        double dDist = calcDist(aadX[i],aadMu[k],nD);

	        if(dDist < dMinDist){
	            nMinK = k;
	            dMinDist = dDist;
	        }
      }

      if(nMinK != anMaxZ[i]){
        int nCurr = anMaxZ[i];
	    nChange++;
	    anW[nCurr]--;
	    anW[nMinK]++;
	    anMaxZ[i] = nMinK;

	    /*check for empty clusters*/
	    if(anW[nCurr] == 0){
	        int nRandI =  gsl_rng_uniform_int (ptGSLRNG, nN);
	        int nKI = 0;
	        /*select at random from non empty clusters*/

	        while(anW[anMaxZ[nRandI]] == 1){
	            nRandI =  gsl_rng_uniform_int (ptGSLRNG, nN);
	        }

	        nKI = anMaxZ[nRandI];
	        anW[nKI]--;
	        anW[nCurr] = 1;
	        anMaxZ[nRandI] = nCurr;
	    }
    }
    }
    //printf("%d %d\n",nIter,nChange);
    nIter++;
    updateMeans(ptCluster, ptData);
  }

  for(i = 0; i < nN; i++){
    for(k = 0; k < nK; k++){
      ptCluster->aadZ[i][k] = 0.0;
    }
    ptCluster->aadZ[i][anMaxZ[i]] = 1.0;
  }

  performMStepMP(ptCluster, ptData);
  return;
}

double eqnA(int nD, gsl_matrix *ptCovarK, gsl_matrix *ptSigmaK, double *adMuK, double *adMK, double dLDetK, double dNuK, double logd2Pi, double dBetaK,double dNK)
{
    int l = 0;
    gsl_matrix *ptRes  = gsl_matrix_alloc(nD,nD);
    gsl_vector *ptDiff = gsl_vector_alloc(nD); 
    double dT1 = 0.0, dT2 = 0.0, dF = 0.0;
    gsl_vector *ptR = gsl_vector_alloc(nD);
    double dD = (double) nD;
    double dRet = 0.0;
    
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0,ptCovarK,ptSigmaK,0.0,ptRes);

    for(l = 0; l < nD; l++){
	    dT1 += gsl_matrix_get(ptRes,l,l);
    }

    for(l = 0; l < nD; l++){
        gsl_vector_set(ptDiff,l,adMuK[l] - adMK[l]);
    }

    gsl_blas_dsymv (CblasLower, 1.0, ptSigmaK, ptDiff, 0.0, ptR);

    gsl_blas_ddot (ptDiff, ptR, &dT2);

    dF = dLDetK - dNuK*(dT1 + dT2) - dD*(logd2Pi + (1.0/dBetaK));

    dRet = 0.5*dNK*dF;

    gsl_matrix_free(ptRes);
    gsl_vector_free(ptDiff);
    gsl_vector_free(ptR);

    return dRet;
}

double eqnB(int nD, gsl_matrix *ptInvW0, gsl_matrix *ptSigmaK, double* adMK, double dBeta0, double d2Pi, double dLDetK, double dBetaK, double dNuK, double dNu0)
{
    int l = 0;
    double dD = (double) nD;
    gsl_matrix *ptRes  = gsl_matrix_alloc(nD,nD);
    gsl_vector *ptDiff = gsl_vector_alloc(nD);
    gsl_vector *ptR = gsl_vector_alloc(nD); 
    double dT1 = 0.0, dT2 = 0.0, dF = 0.0;
    double dRet = 0.0;
    
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0,ptInvW0,ptSigmaK,0.0,ptRes);

    for(l = 0; l < nD; l++){
        dT1 += gsl_matrix_get(ptRes,l,l);
    }

    for(l = 0; l < nD; l++){
	    gsl_vector_set(ptDiff,l,adMK[l]);
    }

    gsl_blas_dsymv (CblasLower, 1.0, ptSigmaK, ptDiff, 0.0, ptR);

    gsl_blas_ddot (ptDiff, ptR, &dT2);

    dF = dD*log(dBeta0/d2Pi) + dLDetK - ((dD*dBeta0)/dBetaK) - dBeta0*dNuK*dT2 - dNuK*dT1;

    dRet = 0.5*(dF + (dNu0 - dD - 1.0)*dLDetK);

    gsl_matrix_free(ptRes);
    gsl_vector_free(ptDiff);
    gsl_vector_free(ptR);
    
    return dRet;
}


double calcVBL_MP(t_Cluster* ptCluster)
{
    int nN = ptCluster->nN;
    int nK = ptCluster->nK, nD = ptCluster->nD;
    double dBishop1 = 0.0, dBishop2 = 0.0, dBishop3 = 0.0, dBishop4 = 0.0, dBishop5 = 0.0; /*Bishop equations 10.71...*/
    double dD = (double) nD;
    double** aadMu = ptCluster->aadMu, **aadM = ptCluster->aadM, **aadZ = ptCluster->aadZ;
    double* adBeta = ptCluster->adBeta, *adNu = ptCluster->adNu, *adLDet = ptCluster->adLDet, *adPi = ptCluster->adPi;
    double adNK[nK], adRet[nK];
    double d2Pi = 2.0*M_PI, logd2Pi = log(d2Pi), dBeta0 = ptCluster->ptVBParams->dBeta0, dNu0 = ptCluster->ptVBParams->dNu0, dRet = 0.0;
    double dK = 0.0;

    for(int k = 0; k < nK; k++){
        adNK[k] = 0.0;
    }

    /*Equation 10.72*/
    for(int i = 0; i < nN; i++){
        for(int k = 0; k < nK; k++){
            adNK[k] += aadZ[i][k];
            if(adPi[k] > 0.0){
	            dBishop2 += aadZ[i][k]*log(adPi[k]);
            }
        }
    }


    for(int k = 0; k < nK; k++){
        if(adNK[k] > 0.0){
            dK++;
        }
    }

    /*Equation 10.71*/
#pragma omp parallel for
    for(int k = 0; k < nK; k++){
        if(adNK[k] > 0.0){
            adRet[k] = eqnA(nD, ptCluster->aptCovar[k], ptCluster->aptSigma[k], aadMu[k],         aadM[k], adLDet[k], adNu[k], logd2Pi, adBeta[k],adNK[k]);
        }
        else{
            adRet[k] = 0.0;
        }
    }
  
    for(int k = 0; k < nK; k++){
        dBishop1 += adRet[k];
    }
  
 
    /*Equation 10.74*/
#pragma omp parallel for
    for(int k = 0; k < nK; k++){
        if(adNK[k] > 0.0){
            adRet[k] = eqnB(nD, ptCluster->ptVBParams->ptInvW0, ptCluster->aptSigma[k], aadM[k], ptCluster->ptVBParams->dBeta0, d2Pi, adLDet[k], adBeta[k], adNu[k],ptCluster->ptVBParams->dNu0);
        }  
    }

    for(int k = 0; k < nK; k++){
        dBishop3 += adRet[k];
    }

    dBishop3 += dK*ptCluster->ptVBParams->dLogWishartB;

    /*Equation 10.75*/
    for(int i = 0; i < nN; i++){
        for(int k = 0; k < nK; k++){
            if(aadZ[i][k] > 0.0){
	            dBishop4 += aadZ[i][k]*log(aadZ[i][k]);
            }
        }
    }

    /*Equation 10.77*/
    for(int k = 0; k < nK; k++){
        if(adNK[k] > 0.0){
            dBishop5 += 0.5*adLDet[k] + 0.5*dD*log(adBeta[k]/d2Pi) - 0.5*dD - dWishartExpectLogDet(ptCluster->aptSigma[k], adNu[k], nD);
        }
    }

    dRet = dBishop1 + dBishop2 + dBishop3 - dBishop4 - dBishop5;

    return dRet;
}

void calcZ_MP(t_Cluster* ptCluster, t_Data *ptData){
    double **aadX = ptData->aadX, **aadZ = ptCluster->aadZ;
    int nK = ptCluster->nK, nD = ptCluster->nD, nN = ptData->nN;
    double dD = (double) nD;
    double** aadM = ptCluster->aadM, *adPi = ptCluster->adPi;

#pragma omp parallel for
    for(int i = 0; i < nN; i++){
        double dMinDist = DBL_MAX;
        double dTotalZ  = 0.0;
        double dNTotalZ = 0.0;
        double adDist[nK];
        int k = 0, l = 0;
        gsl_vector *ptDiff = gsl_vector_alloc(nD);
        gsl_vector *ptRes = gsl_vector_alloc(nD);
    

        for(k = 0; k < nK; k++){
            if(adPi[k] > 0.){
                /*set vector to data point*/
                for(l = 0; l < nD; l++){
                    gsl_vector_set(ptDiff,l,aadX[i][l] - aadM[k][l]);
                }

                gsl_blas_dsymv (CblasLower, 1.0, ptCluster->aptSigma[k], ptDiff, 0.0, ptRes);

                gsl_blas_ddot (ptDiff, ptRes, &adDist[k]);

                adDist[k] *= ptCluster->adNu[k];

                adDist[k] -= ptCluster->adLDet[k];

                adDist[k] += dD/ptCluster->adBeta[k];

                if(adDist[k] < dMinDist){
                    dMinDist = adDist[k];
                }
            }
        }

        for(k = 0; k < nK; k++){
            if(adPi[k] > 0.){
                aadZ[i][k] = adPi[k]*exp(-0.5*(adDist[k]-dMinDist));
                dTotalZ += aadZ[i][k];
            }
            else{
                aadZ[i][k] = 0.0;
            }
        }

        for(k = 0; k < nK; k++){
            double dF = aadZ[i][k] / dTotalZ;
            if(dF < MIN_Z){
                aadZ[i][k] = 0.0;
            }
            dNTotalZ += aadZ[i][k];
        }
        if(dNTotalZ > 0.){
            for(k = 0; k < nK; k++){
                aadZ[i][k] /= dNTotalZ;
            }
        }
        gsl_vector_free(ptRes);
        gsl_vector_free(ptDiff);
    }

    return;
}

void gmmTrainVB_MP(t_Cluster *ptCluster, t_Data *ptData)
{
    int i = 0, k = 0,nIter = 0;
    int nN = ptData->nN, nK = ptCluster->nK;
    /*change in log-likelihood*/
    double dLastVBL = 0.0, dDelta = DBL_MAX;
    double   **aadZ = ptCluster->aadZ;
    int    nMaxIter = ptCluster->nMaxIter;
    double dEpsilon = ptCluster->dEpsilon;
    FILE   *ofp = NULL;

    if(ptCluster->szCOutFile){
	    ofp = fopen(ptCluster->szCOutFile,"w");
	    if(!ofp){
		    fprintf(stderr, "Failed to open file %s in gmmTrainVB\n",ptCluster->szCOutFile);
		    fflush(stderr);
	    }
    }

    /*calculate data likelihood*/
    calcZ_MP(ptCluster,ptData);
    ptCluster->dVBL = calcVBL_MP(ptCluster);

    while(nIter < nMaxIter && dDelta > dEpsilon){

        /*update parameter estimates*/
        performMStepMP(ptCluster, ptData);

        /*calculate responsibilities*/
        calcZ_MP(ptCluster,ptData);

        dLastVBL = ptCluster->dVBL;
        ptCluster->dVBL = calcVBL_MP(ptCluster);
        dDelta = fabs(ptCluster->dVBL - dLastVBL);

        fprintf(stderr,"%d,%f,%f\n",nIter, ptCluster->dVBL, dDelta);
        fflush(stderr);
        if(ofp){
    	    fprintf(ofp,"%d,%f,%f,",nIter, ptCluster->dVBL, dDelta);
    	    for(k = 0; k < nK-1; k++){
    		    fprintf(ofp,"%f,",ptCluster->adPi[k]);
    	    }
    	    fprintf(ofp,"%f\n",ptCluster->adPi[nK - 1]);
	        fflush(ofp);
        }
        nIter++;
    }

    if(ofp){
	    fclose(ofp);
    }

    /*assign to best clusters*/
    for(i = 0; i < nN; i++){
        double dMaxZ = aadZ[i][0];
        int    nMaxK = 0;
        for(k = 1; k < nK; k++){
            if(aadZ[i][k] > dMaxZ){
	            nMaxK = k;
	            dMaxZ = aadZ[i][k];
            }
        }
        ptCluster->anMaxZ[i] = nMaxK;
    }

    return;
}

void calcCovarMatrices(t_Cluster *ptCluster, t_Data *ptData)
{
    int i = 0, j = 0, k = 0, l = 0, m = 0;
    int nN = ptData->nN, nK = ptCluster->nK, nD = ptData->nD;
    double **aadZ = ptCluster->aadZ,**aadX = ptData->aadX;
    double *adPi = ptCluster->adPi, **aadCovar = NULL;
    double dN = (double) nN;

    aadCovar = (double **) malloc(nD*sizeof(double*));
    if(!aadCovar)
        goto memoryError;

    for(i = 0; i < nD; i++){
        aadCovar[i] = (double *) malloc(nD*sizeof(double));
        if(!aadCovar[i])
            goto memoryError;
    }


    for(k = 0; k < nK; k++){ /*loop components*/
        double*     adMu          = ptCluster->aadMu[k];
        gsl_matrix  *ptSigmaMatrix = ptCluster->aptSigma[k];
        /*recompute mixture weights and means*/
        for(j = 0; j < nD; j++){
            adMu[j] = 0.0;
            for(l = 0; l < nD; l++){
	            aadCovar[j][l] = 0.0;
            }
            /*prevents singularities*/
        aadCovar[j][j] = MIN_COVAR;
    }

    /* compute weight associated with component k*/
    adPi[k] = 0.0;
    for(i = 0; i < nN; i++){
        adPi[k] += aadZ[i][k];
        for(j = 0; j < nD; j++){
	        adMu[j] += aadZ[i][k]*aadX[i][j];
        }
    }
    /*normalise means*/
    for(j = 0; j < nD; j++){
        adMu[j] /= adPi[k];
    }

    /*calculate covariance matrices*/
    for(i = 0; i < nN; i++){
        double adDiff[nD];

        for(j = 0; j < nD; j++){
	        adDiff[j] = aadX[i][j] - adMu[j];
        }

        for(l = 0; l < nD; l++){
	        for(m = 0; m <=l ; m++){
	            aadCovar[l][m] += aadZ[i][k]*adDiff[l]*adDiff[m];
	        }
        }
    }

    for(l = 0; l < nD; l++){
        for(m = l + 1; m < nD; m++){
	        aadCovar[l][m] = aadCovar[m][l];
        }
    }

    for(l = 0; l < nD; l++){
        for(m = 0; m < nD; m++){
	        aadCovar[l][m] /= adPi[k];
	        gsl_matrix_set(ptSigmaMatrix, l, m, aadCovar[l][m]);
        }
    }

    adPi[k] /= dN; /*normalise weights*/
  }
  /*free up memory*/
  for(i = 0; i < nD; i++){
    free(aadCovar[i]);
  }

  //gsl_matrix_free(ptSigmaMatrix);
  free(aadCovar);

  return;
 memoryError:
  fprintf(stderr, "Failed allocating memory in performMStep\n");
  fflush(stderr);
  exit(EXIT_FAILURE);
}

/*note assuming you are using inverse W matrix*/
double dLogWishartB(gsl_matrix *ptInvW, int nD, double dNu, int bInv)
{
  int i = 0;
  double dRet = 0.0, dT = 0.0;
  double dLogDet = 0.0, dD = (double) nD;
  gsl_matrix* ptTemp = gsl_matrix_alloc(nD,nD);

  gsl_matrix_memcpy(ptTemp, ptInvW);

  dLogDet = decomposeMatrix(ptTemp, nD);

  if(bInv == TRUE){
    dRet = 0.5*dNu*dLogDet;
  }
  else{
    dRet = -0.5*dNu*dLogDet;
  }

  dT = 0.5*dNu*dD*log(2.0);

  dT += 0.25*dD*(dD - 1.0)*log(M_PI);

  for(i = 0; i < nD; i++){
    dT += gsl_sf_lngamma(0.5*(dNu - (double) i));
  }

  gsl_matrix_free(ptTemp);

  return dRet - dT;
}

double dWishartExpectLogDet(gsl_matrix *ptW, double dNu, int nD)
{
  int i = 0;
  double dRet = 0.0, dLogDet = 0.0, dD = (double) nD;
  gsl_matrix* ptTemp = gsl_matrix_alloc(nD,nD);

  gsl_matrix_memcpy(ptTemp, ptW);

  dLogDet = decomposeMatrix(ptW, nD);

  dRet = dD*log(2.0) + dLogDet;

  for(i = 0; i < nD; i++){
    dRet += gsl_sf_psi(0.5*(dNu - (double) i));
  }

  gsl_matrix_free(ptTemp);

  return dRet;
}

void updateMeans(t_Cluster *ptCluster, t_Data *ptData)
{
  int i = 0, j = 0, k = 0;
  int nN = ptData->nN, nK = ptCluster->nK, nD = ptData->nD;
  int *anMaxZ = ptCluster->anMaxZ;
  int *anW    = ptCluster->anW;
  double **aadX = ptData->aadX, **aadMu = ptCluster->aadMu;

  for(k = 0; k < nK; k++){

    for(j = 0; j < nD; j++){
      aadMu[k][j] = 0.0;
    }
  }

  for(i = 0; i < nN; i++){
    int nZ = anMaxZ[i];

    for(j = 0; j < nD; j++){
      aadMu[nZ][j] += aadX[i][j];
    }
  }

  for(k = 0; k < nK; k++){ /*loop components*/

    /*normalise means*/
    if(anW[k] > 0){
      for(j = 0; j < nD; j++){
	    aadMu[k][j] /= (double) anW[k];
      }
    }
    else{
      for(j = 0; j < nD; j++){
	    aadMu[k][j] = 0.0;
      }
    }
  }

  return;
}

double calcDist(double* adX, double *adMu, int nD)
{
  double dDist = 0.0;
  int i = 0;

  for(i = 0; i < nD; i++){
    double dV = adX[i] - adMu[i];
    dDist += dV*dV;
  }

  return sqrt(dDist);
}
