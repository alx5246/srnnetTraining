import numpy as np
cimport numpy as np
from libc.math cimport exp
import cython


#@cython.boundscheck(False) # turn off bounds-checking for entire function
#@cython.wraparound(False)  # turn off negative index wrapping for entire function
#@cython.nonecheck(False)
def updateCifCoeffsSGDcy(double[:, :] dataMatrix, double[:] spikeArray, double[:, :] filterCoefficients, double learningRate):

    # Define all teh variables
    cdef int numbIts = dataMatrix.shape[0]
    cdef int numbCof = dataMatrix.shape[1]
    cdef double spikeVal
    cdef double predictionOutput

    # Iterate over all the stored data points
    for i in range(numbIts):
        # Generate probabilities
        predictionOutput = 0.0
        for j in range(numbCof):
            predictionOutput += dataMatrix[i,j]*filterCoefficients[j,0]
        predictionOutput = 1./(1. + exp(-1.*predictionOutput))
        spikeVal = learningRate*(spikeArray[-(numbIts-i)] - predictionOutput)
        for j in range(numbCof):
            filterCoefficients[j,0] += spikeVal*dataMatrix[i,j]

    return filterCoefficients
