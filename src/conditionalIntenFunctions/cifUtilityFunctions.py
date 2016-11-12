# A. Lonsberry
# September 2016 - Taken from sandbox (condIntenFunc_0_3.py)
#
# DESCRIPTION
#   Here included are a number of functions that are outside the class objects that are however necessary.


import numpy
import numpy.matlib


########################################################################################################################
# Function Generators
########################################################################################################################
# The CIFs are composed of a set of fucntions weighted by coefficients. We learn the coefficients via SGD or an
# equivelent approach. Here we actually create the numerical versions of these functions which will be convolved with
# a set of spike trains.

def truncatedGaussianFilterGenerator(a, b, c, timeLength, dt, truncation):
    """
    DESCRIPTION
    Create a truncated Gaussian-type filter. Things are truncated so when recorded spike trains are filtered, spikes
    forward in time do not get added to the filter as necessary, and thus keeping data points causel.

    Input 'truncation' is measured in seconds and indicates the first piont in time where the data should be recorded.
    Thus if 'timeLength' is 1.0 seconds which means the Gaussian filter spans from -1 to 1 seconds, and the truncation
    is 0.0, this means all the points in the fitler vector >0 are made to 0.

    IMPORATANT TO NOTE, that because we have a truncation, we are implicitly storing a time-delay in the filter!

    :param a:
    :param b:
    :param c:
    :param timeLength:
    :param dt:
    :param truncation:
    :return:
    """
    # Use numpy.arange and include the last time-step. NOTE, because of rounding errors we have dumped
    # numpy.arange() so not I am using numpy.linspace, but have to be careful not to include the last point! Also we
    # have to be careful using int() as int(6.9999) is 6 and not 7... so round before hand.
    # timeArray = numpy.arange(start=startTime-backupTime, stop=stopTime, step=dt)
    #xvals = numpy.arange(start=-1. * timeLength, stop=timeLength + dt, step=dt)  # Make sure we add the last time-step by adding dt to 'stop'
    numbTimeSteps = int(numpy.round((timeLength*2.0 + dt) / dt)) #We add dt in the summation because we want to include a spot for time=0.
    xvals = numpy.linspace(-1. * timeLength, timeLength, numbTimeSteps)
    convVector = a * numpy.exp(-1. * (((xvals - b) ** 2.) / (2. * c ** 2.)))
    # We want the value at the time = truncation to have a value, but everything greater has to be equal to zero. And
    # more specifically we just get rid of those extra values by deleting them from the array
    convVector = convVector[xvals<=truncation]
    # When later calling numpy.convole, the vector get flipped along the axis, so we pre-flip here so when called it
    # gets flipped back as expected
    convVector = convVector[::-1]
    #print('hey')
    return convVector


########################################################################################################################
# Convolutionl Functions
########################################################################################################################
# In order to find the output of a CIF, we have to convolve some function with a spike train. Below are the function(s)
# required to perform this operation.

def convolveSpikeTimesWthBackUp(convVector, spikeTimes, startTime, stopTime, dt, backupTime=0):
    """
    DESCRIPTION
    Convolve a boolean spike array from either (1) startTime to stopTime or (2) startTime-backupTime to stopTime. In
    either case, the output vector is only going to be over a time interval from [startTime, stopTime]. The inclusion of
    backupTime, is to make the convolution more accurate at the edges.

    The "stopTime" is NOT included, thus the the output is valid over [startTime, stopTime) or [startTime-backupTime,
    stopTime)

    :param convVector: 1D numpy.array, the filter that is convolved with a boolean version of the spike train to produce
        resultant output.
    :param spikeTimes: 1D numpy.array, spike times given in increasing order
    :param startTime: double, the time at which to start tracking, NOTE please pay attention to how rounding to
        time-steps occurs below
    :param stopTime: double, the time at which to stop tracking, NOTE this last time-step is not included
    :param dt: double, the interval of time that is considered a time-step
    :param backupTime: double, the time that we must include backwards. In the case for how we are calculating filters,
        the value for this backupTime is = truncation - (-timeLength) (see self.addFilter() method). Simply this is the
        receding backwards time horizon over which we include spikes from such that the filtered result is correct
        over the time interval we care about.
    :return: [filteredArray, timeArray] where filteredArray is a vector of filtered results at each time-step that
        over the the doamian [startTime, stopTime), and timeArray is also a 1D numpy.array of equal length but
        each part of the vector includes the corresponding time-step value.
    """
    # Use numpy.arange and DO NOT include the last time-step. NOTE, because of rounding errors we have dumped
    # numpy.arange() so not I am using numpy.linspace, but have to be careful not to include the last point! Also we
    # need to be careful with int(), as for example int(5.99999) = 5 and not 6.
    #timeArray = numpy.arange(start=startTime-backupTime, stop=stopTime, step=dt)  # This is going to be for the output
    numbTimeSteps = int(numpy.round((stopTime-dt-startTime)/dt))
    timeArray = numpy.linspace(startTime, stopTime-dt, numbTimeSteps+1) # We have to add 1 to the number of time-steps or else this would not work as intendted (could also remove -dt in line above)
    # Make an array of 0s and 1s to indicate non-spikes and spikes.
    boolArray = numpy.zeros(int(numpy.round((stopTime-(startTime-backupTime))/dt))) # Use round before using int to handle floating point inconsistencies, remember we do not want to inlcude the last time-step
    indsArray = ((spikeTimes - (startTime-backupTime))/dt)
    indsArray = numpy.round(indsArray) # Handle any rounding error before using int() later.
    # Make sure if there is a spike that has a spike-time of end-time, we make it slightly less so it rounds down
    indsArray[indsArray >= boolArray.size] -= 1
    boolArray[indsArray.astype(int)] = 1.
    outputArray = numpy.convolve(boolArray, convVector, mode='valid')
    numbStepsToGrab = int(numpy.round((stopTime-startTime)/dt)) # Use round before using int to handle floating point inconsistencies
    finalOutputArray = outputArray[-numbStepsToGrab:]
    timeArray = timeArray[-numbStepsToGrab:]
    return [finalOutputArray, timeArray]