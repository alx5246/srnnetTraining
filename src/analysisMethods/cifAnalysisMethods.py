# A.Lonsberry
# November 2016
#
# DESCRIPTION
#   Here I am going to have analysis methods specific to Conditional Intensity Functions (CIFs) or more generally
#   comparing or measuring the performance of said CIFs against real network dynamics as measured. Methods here may
#   also be related to quantizing the computational complexity of the CIF and finding equitable trade offs between
#   performance and complexity.

import numpy
import math

def GOFmeth0(cifList, timeArray, dt, spikeTrainList):
    """
    DESCRIPTION
    The goodness of fit method here presented is taken from Haslinger et al. "Discrete Time Rescaling Theorem: ... "
    (2010). This method is used to compare the spike trains and the CIFs that predict them.
    :param cifList: list of 1D numpy.array() with values of lamba, where the bin of the array has time-step of dt
    :param timeArray: 1D numpy.array() with the values of time (in units of seconds) at which the bin begins (ie
           timeArray[0] = 0.0 seconds means the first bin of time is recorded from time t given 0.0<=t to time<0.0+dt,
           that is the end of the bin is not included)
    :param dt: time-step that the CIFs functions are sampled at (in units of seconds).
    :param spikeTrainList: list of 1D numpy.array() that indicates the spike times.
    :return:
    """

    # Iterate over list of CIFs (1D numpy.arrays all stored in cifList)
    for i,cif in enumerate(cifList):

        # Step 1, find which bins (indices) the spikes belong to. We will store this as a new array of the same size as
        # the array that holds the spike times.
        spikeTimes = spikeTrainList[i]
        if spikeTimes.shape[0]>0: # Make sure their are spike times to look at!
             spikeInds = numpy.floor((spikeTimes - timeArray[0]) / dt) # The bin of time which the spike occurs.

        # Step 2, generate a new time series of discrete q(k) values, where k is the index and q(k) = -log(1-p(k)) where
        # p(k) is the probability at the time-step and is found as p(k)=~dt*lamda(k), where lambda are the values stored
        # in the numpy.arrays given in the cifList.
        q = -1.*numpy.log(1.-cif*dt)

        # Step 3, for each interspike-interval j, we calculate the rescaled ISI deemed as zeta(j)
