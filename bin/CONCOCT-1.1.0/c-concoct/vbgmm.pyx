"""
vbgmm.pyx

simple cython wrapper for variational Gaussian mixture model in C 

"""

import cython

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern void c_vbgmm_fit (double* adX, int nN, int nD, int nK, int seed, int* anAssign, int nThreads, int nIter)
@cython.boundscheck(False)
@cython.wraparound(False)

def fit(np.ndarray[double, ndim=2, mode="c"] xarray not None, nClusters, seed, threads, piter):
    """
    fit (xarray, nClusters, seed, threads)

    Takes a numpy array xarray as input, fits the vbgmm using nClsuters initial clusters

    param: xarray -- a 2-d numpy array of np.float64
    param: nClusters -- an int, number of start clusters
    param: seed -- an int, the random seed
    param: threads -- int, the number of threads to use
    param: piter -- int, the number of VB iterations to use
    """
    cdef int nN, nD, nK, nThreads, nIter
        
    nN, nD = xarray.shape[0], xarray.shape[1]

    nK = nClusters

    nIter = piter

    nThreads = threads

    cdef np.ndarray[int, ndim=1,mode="c"] assign = np.zeros((nN), dtype=np.intc)
    
    c_vbgmm_fit (&xarray[0,0], nN, nD, nK, seed, &assign[0], nThreads, nIter)

    return assign
