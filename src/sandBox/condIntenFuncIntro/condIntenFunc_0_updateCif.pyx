import numpy as np
cimport numpy as np
from libc.math cimport exp
import cython

DTYPE = np.float
DTYPEint = np.int

ctypedef np.float_t DTYPE_t
ctypedef np.int_t DTYPEint_t

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.nonecheck(False)
def updateCifCoeffsSGDcy(np.ndarray[DTYPE_t, ndim=2] dataMatrix, np.ndarray[DTYPE_t, ndim=1] spikeArray, np.ndarray[DTYPE_t, ndim=2] filterCoefficients, float learningRate):

    assert dataMatrix.dtype == DTYPE and filterCoefficients.dtype == DTYPE and spikeArray.dtype == DTYPE

    # Define all teh variables
    cdef int numbIts = dataMatrix.shape[0]
    cdef int numbCof = dataMatrix.shape[1]
    cdef DTYPE_t predictionOutput

    cdef DTYPE_t spikeVal
    cdef DTYPE_t value

    # Iterate over all the stored data points
    for i in range(numbIts):
        # Generate probabilities
        predictionOutput = 0.0
        for j in range(numbCof):
            predictionOutput += dataMatrix[i,j]*filterCoefficients[j,0]
        predictionOutput = 1./(1. + exp(-1.*predictionOutput))
        # spikeVal = spikeArray[-(numbIts-i)]
        spikeVal = learningRate*(spikeArray[-(numbIts-i)] - predictionOutput)
        for j in range(numbCof):
            # value = filterCoefficients[j,0] + learningRate*(spikeVal - predictionOutput)*dataMatrix[i,j]
            # value = filterCoefficients[j,0] + spikeVal*dataMatrix[i,j]
            # filterCoefficients[j,0] = value
            filterCoefficients[j,0] = filterCoefficients[j,0] + spikeVal*dataMatrix[i,j]

    return filterCoefficients
