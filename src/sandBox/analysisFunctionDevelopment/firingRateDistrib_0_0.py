# A. Lonsberry
# May 2016
# (Using python 3.4.4, and brian2 2.0rc)
#
# DESCRIPTION
#   Here develop a function/method to produce a smooth firing-rate approximation of neuron spike output.
#
#   Here we demonstrate a method for looking at the firing-rate of a single neuron vs time. In particular, the method
#   here will take the spike-train (that could be recorded for example by using a 'SpikeMonitor' from brian2), and
#   output a smooth firing-rate approximation.
#
#   Here the particular methods are found, in another file (firingRateDistrib_0_1.py) we will actually be calling and
#   testing the methods.

import numpy

def instant_firing_rate(spikeTrain, startTime, endTime, dt=.001, filterType = 'Gaussian', filterLength=.5, var = .25):
    """
    FUNCTION DESCRIPTION
        This function takes as input a 1-D numpy.array of spike times, and outputs a firing frequency vector for each
        time-step.

        This is an example function, and thus has some limitaions. For example it can only handle 1D input. A more
        complex version of this function would likely be able to handle N-D input (spike trains from multiple neurons)
        simultaneously.

    :param spikeTrain: 1D numpy.array, units are seconds, neuron spike times stored in an numpy.array
    :param startTime: scalar, units are seconds, the time at which we want to start approximating the firing-rate
    :param endTime: scalar, units are seconds, the time at which we want to stop approximating the firing-rate
    :param dt: scalar, units are seconds, dt=0.0001 is the standard
    :param filterType: string, the type of smoothing filter used, options are either 'Gaussian' or 'Square'
    :param filterLength: scalar, units are seconds, the lenght of the smoothing filter
    :param var: scalar, units : 1, the variance of the Gaussian function used if Guassian filter is applied
    :return: 1X2 list, [convolvedVector, timeVector], both convolvedVector and timeVector are of numpy.arrray types
    """

    #-------------------------------------------------------------------------------------------------------------------
    # Create Convolution Vectors
    #-------------------------------------------------------------------------------------------------------------------_
    # The mechanism to find the firing rate of the spike train, is to convolve a vector of ones (spikes) and zeros
    # (no-spikes) with a filter. In the case here the filter will be either square or Gaussin.

    # Below we create the vectors of value we want to convolve against a binary spike train (1 where spike, 0 where no
    # spike). We have different types of filters to do this, the standard is the "Gaussian" type.
    if filterType=='Gaussian':

        # Make an array we can give the Gaussain function
        timeArray = numpy.linspace(start=0-filterLength/2, stop=filterLength/2, num=int(filterLength/dt), endpoint=True)
        # Now make the Gaussian convolution vector
        filterVector = numpy.exp(-(timeArray**2)/(2*(var**2)))

    else :

        # Make an array of equal values
        filterVector = numpy.ones(int(filterLength / dt))

    # Find the normalizing coefficient, and normalize the firing frequency.
    normCoef = (1./dt)*(1./numpy.sum(filterVector))
    filterVector *= normCoef

    #-------------------------------------------------------------------------------------------------------------------
    # Create Boolean Spike Train
    #-------------------------------------------------------------------------------------------------------------------
    # We cannot convolve the array of spike-times directly, rather we create an array of ones (spikes) and zeros
    # (no-spikes)

    # Make spike train into somthing we can convolve with
    numberTimeSteps   = int((endTime-startTime)/dt) #New array of zeros and ones to represent the spike train
    spikeTrainBoolean = numpy.zeros(numberTimeSteps)

    # Grab the spike times we care about
    #condition = spikeTrain>=startTime, spikeTrain<endTime
    #print(len(condition[0]))
    #print(condition[1])
    #print(spikeTrain.shape)
    #spikeTrainSliced = numpy.extract(condition, spikeTrain)
    spikeTrainSliced = spikeTrain[(spikeTrain>startTime) & (spikeTrain<endTime)]
    print(spikeTrainSliced)
    spikeTrainMapped  = spikeTrainSliced/dt # Here we cast the values
    print(spikeTrainMapped)

    # Add spikes to the boolean array
    for aSpike in spikeTrainMapped:
        index = int(aSpike - startTime/dt)
        print(index)
        spikeTrainBoolean[index] += 1

    #-------------------------------------------------------------------------------------------------------------------
    # Convolve
    #-------------------------------------------------------------------------------------------------------------------

    #Convolve the boolean vector with the filter we would like to use.
    convolvedVector = numpy.convolve(spikeTrainBoolean, filterVector, mode='same')

    #-------------------------------------------------------------------------------------------------------------------
    # Time Vector Output
    #-------------------------------------------------------------------------------------------------------------------
    # I need to make an output vector for the time-axis
    timeVector = numpy.linspace(startTime,endTime,num=convolvedVector.shape[0])

    #-------------------------------------------------------------------------------------------------------------------
    # Output
    #-------------------------------------------------------------------------------------------------------------------
    return [convolvedVector, timeVector]

