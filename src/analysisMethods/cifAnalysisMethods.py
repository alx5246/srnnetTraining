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

# I had a difficult time getting things to import, so now I have however, go about fixing this by going up a directoy
# which I think I can do in the following way.
import os
import sys
sys.path.append(os.path.abspath("../"))

def GOFmeth0(cifList, timeArray, dt, spikeTrainList):
    """
    DESCRIPTION
    The goodness of fit method here presented is taken from Haslinger et al. "Discrete Time Rescaling Theorem: ... "
    (2010). This method is used to compare the spike trains and the CIFs that predict them. Given an input list of
    CIFs, and the associated list of spikeTrains that may or may not match the CIFs, we where generate the information
    to compare both CIF to spike-train. The output of this function is a list of numpy.arrays, a numpy.array for each
    CIF (and corresponding spike-train) given. These numpy.arrays are filled with sorted (ascending) values of
    rescaled interspike intervals. If the values are uniformly distributed (perhaps use a histogram to look at and a KS
    test to measure), then the spike-trains are likely from the CIF.
    :param cifList: list of 1D numpy.array() with values of lamba, where the bin of the array has time-step of dt
    :param timeArray: 1D numpy.array() with the values of time (in units of seconds) at which the bin begins (ie
           timeArray[0] = 0.0 seconds means the first bin of time is recorded from time t given 0.0<=t to time<0.0+dt,
           that is the end of the bin is not included)
    :param dt: time-step that the CIFs functions are sampled at (in units of seconds).
    :param spikeTrainList: list of 1D numpy.array() that indicates the spike times.
    :return: lsit of numpy.arrays, each numpy array will have the transformed y_i and stored (ascending) interspike intervals.
    """

    outputList = []

    # Iterate over list of CIFs (1D numpy.arrays all stored in cifList)
    for i,cif in enumerate(cifList):

        # Step 1, find which bins (indices) the spikes belong to. We will store this as a new array of the same size as
        # the array that holds the spike times.
        spikeTimes = spikeTrainList[i]
        if spikeTimes.shape[0]>0: # Make sure their are spike times to look at!
            spikeIndsDouble = numpy.floor((spikeTimes - timeArray[0]) / dt) # The bin of time which the spike occurs.
            spikeInds = spikeIndsDouble.astype(int)
            # Step 2, generate a new time series of discrete q(k) values, where k is the index and q(k) = -log(1-p(k)) where
            # p(k) is the probability at the time-step and is found as p(k)=~dt*lamda(k), where lambda are the values stored
            # in the numpy.arrays given in the cifList.
            q = -1.*numpy.log(1.-cif*dt)

            # Step 3, for each interspike-interval j, we calculate the rescaled ISI deemed as zeta(j)
            # To do this we must generate a sequence of random values called 'delta' for each j interval, so we need to
            # know how many inter-spike intervals there are
            nSpikeIntervals = spikeInds.shape[0] - 1 # The reason we subtract 1 is that we care about the interval, and there are 1 less intervals than spikes
            deltaRand = numpy.random.random((nSpikeIntervals,1))
            # Now we calculate a simplified version of the delta values (some of the equations appear to cancel out)
            expComp = 1.-(numpy.exp(-q[spikeInds[1:spikeInds.size]]))
            delta = -1. * numpy.log(1-deltaRand*(expComp))

            zeta = numpy.zeros(spikeInds.shape[0] - 1)
            for j, deltaj in enumerate(delta):
                dummyA = q[spikeInds[j]+1:spikeInds[j+1]]
                zeta[j] = numpy.sum(q[spikeInds[j]+1:spikeInds[j+1]]) + deltaj[0]

            # Step 4, make final transform to yi
            y = 1 - numpy.exp(-zeta)

            # Now if the statistical model is accurate the y vlaues will be uniformly distributed, and therefor the y
            # can be used to make a KS or differential KS plot. In order to find the plot, we order the values from
            # smallest to largest and then compare against a uniform CFD.
            ySort = numpy.sort(y)
            outputList.append(ySort)

        else:
            # If there are not interspike intervals... send out an empty array
            outputList.append(numpy.array([]))

        return outputList


def testGOFmeth0():
    """
    DESCRIPTION
    This is simply a method to test the method I formally created.
    :return:
    """
    # Import the module we require in order to make som arbitrary spike trains.
    from spikeDataGenerators import analogToSpikes
    import matplotlib.pyplot as plt
    import scipy.stats
    # Generate some random data that will serve as a surrogate CIF (lambda values). This data is hardcoded
    xVals = numpy.random.random((10000000, 1)) * 1.
    timeValues = numpy.linspace(0., 9999.999, 10000000)
    # Generate a Point-Process using Time-Rescaling, that is we use our other code to generate spikes from our randomly
    # generated analog data.
    spikePntPrcs = analogToSpikes.genSpikesWithTimeRescaling(xVals, .001, 1.)
    spikePntPrcs = numpy.asarray(spikePntPrcs[0])  # Grab the first numpy.array in the list
    # Plot the spikes and xValues
    spikePntPrcsPlt = numpy.zeros(spikePntPrcs.shape[0] * 3)
    spikePntPrcsPltTimes = numpy.zeros(spikePntPrcs.shape[0] * 3)
    index = 0
    for i in range(spikePntPrcs.shape[0]):
        spikePntPrcsPlt[index] = 0.
        spikePntPrcsPltTimes[index] = spikePntPrcs[i]
        index += 1
        spikePntPrcsPlt[index] = 1.
        spikePntPrcsPltTimes[index] = spikePntPrcs[i]
        index += 1
        spikePntPrcsPlt[index] = 0.
        spikePntPrcsPltTimes[index] = spikePntPrcs[i]
        index += 1
    plt.figure(1)
    plt.plot(spikePntPrcsPltTimes,spikePntPrcsPlt)
    # plt.show()
    # Pass to our method we want to test to see if it works. We will know if it works by seeing if the values that come
    # out are drwn from a uniform distribution. That means the probability of getting of one value is the same as the
    # probabiltiy of getting another value over the domain.
    yTrnsList = GOFmeth0([xVals], timeValues, dt=.001, spikeTrainList=[spikePntPrcs])
    plt.figure(2)
    plt.hist(yTrnsList[0], bins=50)
    # Now using a KS test, see if the output is like a uniform distribution, the more like a uniform distribution, the
    # more the spikes are likely to be generated from the time-varying input.
    print(scipy.stats.kstest(yTrnsList[0], 'uniform'))
    plt.show()


########################################################################################################################
# Internal Unit Testing
########################################################################################################################
# I want to do some very simple unit testing on the functions above.

if __name__ == "__main__":

    # TEST GOFmeth0()
    testGOFmeth0()



