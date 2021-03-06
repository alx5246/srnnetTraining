# A. Lonsberry
# June 2017
#
# DESCRIPTION
#   This python file will have method(s) to create arrays of spike times that can be used as input to some simulation.

import math
import random
import numpy

def poissonSpikeGen(rateArray, timeArray):
    '''
    DESCRIPTION
    Creates a 1-D python list of spike times between timeArray[0] and timeArray[-1]. The times are created by Poisson
    process given a rate function(s) that are input. In the case here the rate functions are specified for
    time-intervals given in the input 'timeArray'.

    EXAMPLE
    Given (1) three time intervals 1.0 to 2.5 seconds, 2.5 seconds to 3.5 seconds, and 3.5 seconds to 7.0 seconds, and
    (2) the rates over those time-intervals 1.0 Hz, 3.0Hz, and 4.0Hz, the input to this functoin is given as
    rateArray = [1., 3., 4.] and timeArray = [1.0, 2.5, 3.5, 7.0].

    :param rateArray: python 1-D list, of rate values over time intervals
    :param timeArray: pthon 1-D list, with time intervals of constant rates
    :return: python 1-D list, of spike times
    '''

    #Make initial element, which will at the end be removed! This is needed however to create the rest of the list.
    spikeTimes = [timeArray[0]]

    # Make our first sample
    countIndex    = 0
    tarVal        = -1. * math.log(1. - random.random())
    if tarVal <= .1:
        tarVal = .11

    intermedTime  = 0.
    intermedScore = 0.

    # Iterative until having gone past time limit
    while countIndex < len(rateArray):

        # Find delta-time
        delTime = (tarVal - intermedScore) / rateArray[countIndex]

        if (spikeTimes[-1] + delTime + intermedTime) > timeArray[countIndex+1]:
            # We do not crate a spike, we need to go into the next rate interval and accumulate more time.

            dummyTime     =  min( timeArray[countIndex+1]-spikeTimes[-1] , timeArray[countIndex+1]-timeArray[countIndex] )
            intermedTime  = intermedTime + dummyTime
            intermedScore = intermedScore + rateArray[countIndex]*dummyTime
            countIndex += 1

        else:
            # We need to create a new spike and draw a new values
            spikeTimes.append(spikeTimes[-1] + delTime + intermedTime)
            tarVal = -1. * math.log(1. - random.random())
            if tarVal <= .1:
                tarVal = .11
            intermedTime  = 0.
            intermedScore = 0.

    # Drop the first value
    del spikeTimes[0]

    return spikeTimes


def multiPoissonStreams(numPoiStreams, startTime, stopTime, rateBounds, intervalBounds):
    '''
    DESCRIPTION
    Make a number of input streams, where for example, each input stream could correspond to spikes from some input
    neurons.


    :param numPoiStreams: integer > 0, number of Poisson spike time streams
    :param startTime: double > 0., the smallest spike time allowed
    :param stopTime: double > startTime, the largest spike time allowed
    :param rateBounds: 1X2 list, [minrate, maxrate] assumed to be in Hz
    :param intervalBounds: 1X2 list. [minInterval, maxInterval] assumed to be in Seconds

    :return: 2D numpy.array, number of rows equal total number of spikes, there are two columns. In the first column is
    the index of the neuron that spiked, and in the second column is the time of said spike. The array is sorted by
    time in asscending order.
    '''

    streamIndex = numpy.array([])
    spikeTimes = numpy.array([])

    for i in range(numPoiStreams):

        #Create lists of time intervals, and list of firing rates.

        timeIntervals = [startTime]
        rateIntervals = []

        timeSample = random.random() * (intervalBounds[1]-intervalBounds[0]) + intervalBounds[0]
        rateSample = random.random() * (rateBounds[1]-rateBounds[0]) + rateBounds[0]

        while timeIntervals[-1] + timeSample < stopTime:

            newTime = timeIntervals[-1] + timeSample
            timeIntervals.append(newTime)
            rateIntervals.append(rateSample)
            timeSample = random.random() * (intervalBounds[1] - intervalBounds[0]) + intervalBounds[0]
            rateSample = random.random() * (rateBounds[1] - rateBounds[0]) + rateBounds[0]

        #Now make last sample
        timeIntervals.append(stopTime)
        rateSample = random.random() * (rateBounds[1] - rateBounds[0]) + rateBounds[0]
        rateIntervals.append(rateSample)

        #Now give to the spike maker function
        newSpikeTimes = poissonSpikeGen(rateIntervals, timeIntervals)

        #Now add to our ever increasing lists of spikes and spike generators
        streamIndex = numpy.concatenate((streamIndex, numpy.ones(numpy.asarray(newSpikeTimes).size, dtype=int)*i))
        spikeTimes  = numpy.concatenate((spikeTimes, numpy.asarray(newSpikeTimes)))

    #Now with our two numpy arrays filled we want to stick them together and then sort them
    outputArray = numpy.transpose(numpy.vstack((streamIndex, spikeTimes)))
    outputArray = outputArray[outputArray[:, 1].argsort()]

    return outputArray


########################################################################################################################
# Internal Unit Testing
########################################################################################################################

# I want to do some very simple unit testing on the functions above.

if __name__ == "__main__":

    #TEST THE FUNCTIONS

    #Test (1)
    spikeArray = poissonSpikeGen([1.,2.,1.], [1., 3., 4., 6. ] )
    print(spikeArray)
    print('\n')

    #Test (2)
    spikesArray = multiPoissonStreams(3, startTime=0.0, stopTime=10.0, rateBounds=[1., 5.], intervalBounds=[1., 3.])
    print(spikesArray)







