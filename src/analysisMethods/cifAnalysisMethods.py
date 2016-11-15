# A.Lonsberry
# November 2016
#
# DESCRIPTION
#   Here I am going to have analysis methods specific to Conditional Intensity Functions (CIFs) or more generally
#   comparing or measuring the performance of said CIFs against real network dynamics as measured. Methods here may
#   also be related to quantizing the computational complexity of the CIF and finding equitable trade offs between
#   performance and complexity.

import numpy
import numpy.random

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
    :return: Kolmogorov-Smirnov Statistic Dn (see wikipedia... lawl)
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
            # To do this we must generate a sequence of random values called 'delta' for each j interval, so we need to
            # know how many inter-spike intervals there are
            nSpikeIntervals = spikeInds.shape[0] - 1 # The reason we subtract 1 is that we care about the interval, and there are 1 less intervals than spikes
            deltaRand = numpy.random.random(nSpikeIntervals)
            # Now we calculate a simplified version of the delta values (some of the equations appear to cancel out)
            delta = -1. * numpy.log(1-deltaRand*(1-numpy.exp(-q[spikeInds[1:spikeInds.size]])))
            zeta = numpy.zeros(spikeInds.shape[0] - 1)
            for j, deltaj in enumerate(delta):
                zeta[j] = numpy.cumsum(q[spikeInds[j]+1:spikeInds[j+1]]) + deltaj

            # Step 4, make final transform to yi
            y = 1 - numpy.exp(-zeta)

            # Noow if the statistical model is accurate the y vlaues will be uniformly distributed, and therefor the y
            # can be used to make a KS or differential KS plot. In order to find the plot, we order the values from
            # smallest to largest and then compare against a uniform CFD.
            ySort = numpy.sort(y)

