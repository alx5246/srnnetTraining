import numpy as np
cimport numpy as np
from libc.math cimport round
import cython

DTYPE = np.int
ctypedef np.int_t DTYPE_t
DTYPEd = np.double
ctypedef np.double_t DTYPEd_t



#@cython.boundscheck(False) # turn off bounds-checking for entire function
#@cython.wraparound(False)  # turn off negative index wrapping for entire function
#@cython.nonecheck(False)
def convSpikeTimesWithBackup(double[:] convVector, double[:] spikeTimes, double startTime, double stopTime, double dt, double backupTime):

    # Define all the variables
    cdef int numbSpikes = spikeTimes.shape[0]
    cdef int numbTimeSteps = np.int(round((stopTime-startTime)/dt))
    cdef int numbTimeStepsWback = np.int(round(((stopTime-(startTime-backupTime))/dt)))
    cdef int convSize = convVector.shape[0]
    cdef int ind
    cdef int counter = 0
    cdef DTYPEd_t value
    cdef np.ndarray timeArray = np.zeros(numbTimeSteps, dtype=DTYPEd)
    cdef np.ndarray boolArray = np.zeros(numbTimeStepsWback, dtype=DTYPEd)
    cdef np.ndarray outputArray = np.zeros(numbTimeSteps, dtype=DTYPEd)


    # Fill in time values
    for i in range(numbTimeSteps):
        timeArray[i] = startTime + dt*i

    # Fill in spikes
    for i in range(numbSpikes):
        ind = np.int(round((spikeTimes[i] - (startTime-backupTime))/dt))
        if ind>=numbTimeStepsWback:
            ind -= 1
        boolArray[ind] = 1.

    # Perform convolution (remeber to flip the vector backwards)
    for i in range(numbTimeSteps):
        value = 0.0
        for j in range(convVector):
            value += convVector[-j]*boolArray[counter+i]
        outputArray[i] = value

    return [outputArray, timeArray]



