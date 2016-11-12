# A. Lonsberry
# September 2016
#
# DESCRIPTION
#   Methods and/or classes necessary to take analog data and decompose into different components.
#
#   For example, suppose we have analog data along 1 dimension. This may not be well suited as input into a network
#   of neurons. Thus we decompose that axis into seperate set of verlapping regions wherein a neuron fires when the
#   analog input is within this region.
#
#   Included are methods to find coefficients in Gaussians to make sure we have an optimal overlap over some specific
#   domain. Included are methods to decompose 1D data via Gaussians.



import numpy
import numpy.random
import math
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection as pc
from mpl_toolkits.mplot3d import Axes3D




def sumOverlapGaussians(xStart, xStop, gaussianPoints, b, dt=.001):
    """
    DESCRIPTION
    Simply outputs the summation of over-laping Gaussian functions at each point from xStart to xStop. This was mostly
    generated as a method to help check the output of other methods in this .py file.
    :param xStart:
    :param xStop:
    :param gaussianPoints:
    :param b:
    :param dt:
    :return:
    """
    time = numpy.arange(xStart, xStop, dt)
    cfd = numpy.zeros(numpy.shape(time)[0])
    for i in range(time.shape[0]):
        for j in range(numpy.shape(gaussianPoints)[0]):
            cfd[i] = cfd[i] + (1./((b**.5)*math.pi)) * math.exp(-1.*((time[i]-gaussianPoints[j])**2.0)/b)
    return cfd


def minVarOverlappingGaussian(xStart, xStop, nGaus, b, dt=.001):
    """
    DESCRIPTION
        Running a bisection method to find 'b' such that we get non varying cfd between gaussian peaks, meaning we have
        a monotonic function all the way to the center. This ends up working very well and very fast.
    :param xStart: Position of first Gaussian
    :param xStop: Position of second Gaussian
    :param nGaus: Number of Gaussians
    :param b: b starting point
    :param dt: discretizationa long x axis (shouldn't this be dx?)
    :return: [gCenters, bMid] where gCenters is a 1D numpy array that indicates the positions of the Gaussian functions,
             and bMid=2*sigma^2 is the 'optimal' value such that the if we look at the summation of Gaussians at each
             point in the domain, we end up with a function that is just monotonically increasing for the first half of
             the domain, and is monotically decreasing for the second half of the domain.
    """

    # Firstly we are assuming how we are spacing the Gaussians across the domain. Find the Gaussain center values
    # and a time vector as well.
    dStep = (xStop - xStart) / (nGaus - 1)
    gCenters = numpy.arange(xStart, xStop + dStep, dStep)
    time = numpy.arange(xStart, xStop+dt, dt)

    # Find half steps number
    hS = int((numpy.shape(time)[0])/2.)

    # Create the CFD values first
    cfdHigh = numpy.zeros(numpy.shape(time)[0]) # Too smooth
    bHigh = b
    cfdLow = numpy.zeros(numpy.shape(time)[0])  # Too smooth
    bLow = b

    # First find an acceptable value of b for the high values
    isGood = False
    while isGood==False:
        strOut = "bHigh = " + str(bHigh)
        print(strOut)
        # Calculate the highCFD and see if we get undulation or ... not
        cfdHigh = numpy.zeros(numpy.shape(time)[0])  # Too smoot
        for i in range(time.shape[0]):
            for j in range(nGaus):
                cfdHigh[i] = cfdHigh[i] + (1. / ((bHigh * math.pi ** .5) )) * math.exp(-1. * ((time[i] - gCenters[j]) ** 2.0) / bHigh)
        dummyVec = cfdHigh[1:hS] - cfdHigh[0:hS-1]
        # Check if all values are positive!
        if not numpy.all(dummyVec>=0.0):
            bHigh = bHigh*1.05
        else:
            isGood = True

    # Second find an acceptable value of b for the low values
    isGood = False
    while isGood == False:
        strOut = "bLow = " + str(bLow)
        print(strOut)
        # Calculate the highCFD and see if we get undulation or ... not
        cfdLow = numpy.zeros(numpy.shape(time)[0])  # Too smooth
        for i in range(time.shape[0]):
            for j in range(nGaus):
                cfdLow[i] = cfdLow[i] + (1. / ((bLow * math.pi ** .5) )) * math.exp(-1. * ((time[i] - gCenters[j]) ** 2.0) / bLow)
        dummyVec = cfdLow[1:hS] - cfdLow[0:hS - 1]
        # Check if all values are positive!
        if numpy.all(dummyVec >= 0.0):
            bLow = bLow * 0.95
        else:
            isGood = True


    # Now we need to make an intermiediate value!
    cfdMid = numpy.zeros(numpy.shape(time)[0])  # Too smooth
    bMid = (bHigh-bLow)/2.0 + bLow

    #Iterate
    while bHigh-bLow > .0000001:
        strOut = "bHigh-bLow =  " + str(bHigh-bLow )
        print(strOut)
        # Calculate the midCFD and see if we get undulation or ... not
        cfdMid = numpy.zeros(numpy.shape(time)[0])  # Too smooth
        for i in range(time.shape[0]):
            for j in range(nGaus):
                cfdMid[i] = cfdMid[i] + (1. / ((bMid * math.pi ** .5) )) * math.exp(-1. * ((time[i] - gCenters[j]) ** 2.0) / bMid)
        dummyVec = cfdMid[1:hS] - cfdMid[0:hS - 1]
        # Check which one we should replace
        # Check if all values are positive!
        if numpy.all(dummyVec >= 0.0):
            bHigh = bMid
            cfdHigh = cfdMid
        else:
            bLow = bMid
            cfdLow = cfdMid
        # Recalculate b values
        bMid = (bHigh - bLow) / 2.0 + bLow

    return [gCenters, bMid]


def seperateInputValuesAlongGaussians(xvalues, gcenters, b, normalize=False):
    """
    DESCRIPTION
    Takes in some 1D values (use must know the domain they are over), some Gaussians defined by 'gcenters' and 'b', and
    breaks up the input into overlapping parts. Note we define the Guassian in this case as a normal function of the
    form f(x) = 1/((b*pi)**.5) * exp( - ((x-mu)**2)/(b) )
    :param xvalues: 1D numpy.array over the domain we want to seperate out
    :param gcenters: The centers along the 1D domain
    :param b: Part of the Normal function
    :param normalize: Not yet used
    :return: 2D numpy.array [xvalues.size, gcenters.size], teh 1D set of values decomposed into seperate components.
    """
    # Make empty array to fill
    segmentedValues = numpy.zeros([xvalues.size,gcenters.size])
    for i, x in enumerate(xvalues):
        for j, gc in enumerate(gcenters):
            segmentedValues[i, j] = (1./((b*math.pi)**.5)) * math.exp( -((x - gc)**2.0)/b )

    return segmentedValues

########################################################################################################################
# UNIT TESTING
########################################################################################################################
# Everything seems to be working correctly


if __name__ == "__main__":

    import analogToSpikes

    # TRY OUT BISECTION ALGORITHM TO FIND APPROPRIATE GAUSSIANS
    # Use bisection alogithm above to find parameters of Guassian so to create 'optimal' set of overlapping functions.
    #plt.figure(1)
    #xStart = 0.0
    #xStop = 10.
    #cfds = minVarOverlappingGaussian(xStart, xStop, 11, .05, .001)
    #gCenters = cfds[0]
    #b = cfds[1]
    #cfd = sumOverlapGaussians(xStart, xStop, gCenters, b, .001)
    #plt.plot(cfd)
    #print(b)
    #plt.show()

    # TRY OUT THE DECOMPOSITION OF 1D DATA
    # Make some gaussins placed over the 1d space (faster then running bisection above)
    gCenters = numpy.array([0.0, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.])
    b = 1.775
    # Now make some data and send it into the thing to be parsed
    dArray = numpy.linspace(0.0, 5.0, num=5001)
    x = (dArray - 2.5)**2.
    #Now break up parsed data into the correct components
    segmentedValues = seperateInputValuesAlongGaussians(x, gCenters, b)

    #Now lest plot shit
    fig = plt.figure(2)
    ax = Axes3D(fig) # Because I am using older version
    verts = []
    for i in range(gCenters.size):
        segmentedValues[0,i] = 0
        segmentedValues[-1,i] = 0
       # print(list(zip(segmentedValues[:,i],dArray)))
        verts.append(list(zip(segmentedValues[:,i],dArray)))
    poly = pc(verts)
    ax.add_collection3d(poly, gCenters, zdir='y')
    ax.set_xlim3d(0,1)
    ax.set_zlim3d(0,5)
    ax.set_ylim3d(0,10)
    #plt.show()

    # TRY TAKING ANALOG TO SPIKES (Based on Poisson ReScaling bit) (Or based on simple linear scaling)
    # spikeTrains = analogToSpikes.genSpikesWithTimeRescaling(segmentedValues, dt=.001, scaling=10.)
    spikeTrains = analogToSpikes.genSpikesLinearly(segmentedValues, dt=.001, scaling=20.)
    # Now lets try to plot this bitch, iterate over the different spike trains
    verts = []
    for i in range(len(spikeTrains)):
        if len(spikeTrains[i])>0:
            # Convert the list to a numpy array so we can do stuff there
            aSpikeTrain = numpy.asarray(spikeTrains[i])
            #
            x = numpy.zeros(aSpikeTrain.size*3)
            y = numpy.zeros(aSpikeTrain.size*3)
            ind = 0
            for j in range(aSpikeTrain.size):
                x[ind] = aSpikeTrain[j]
                y[ind] = 0.
                ind += 1
                x[ind] = aSpikeTrain[j]
                y[ind] = 1.
                ind += 1
                x[ind] = aSpikeTrain[j]
                y[ind] = 0.
                ind += 1
            #
            print(list(zip(x,y)))
            verts.append(list(zip(y,x)))
        else:
            # Need to give something to those without spikes or graph does not work out
            x = numpy.array([0.,.001])
            y = numpy.array([.001,.001])
            verts.append(list(zip(y, x)))
    poly = pc(verts)
    ax.add_collection3d(poly, gCenters, zdir='y')
    plt.show()






