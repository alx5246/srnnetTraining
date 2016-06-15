# A. Lonsberry
# June 2017
#
# DESCRIPTION
#   This python file will have methods to create and train the conditional intensity functions that we use to predict
#   the spiking of the neurons.
#
#   In the simplest version we will have a single filter, we look at its value for different delays, this means we only
#   have to run a convolution operation for each neuron once!!! This should really speed up things.

import numpy
import matplotlib.pyplot as plt

########################################################################################################################
# Create Individual Elements of CIF
########################################################################################################################


def gaussianFilterGenerator(a, b, c, dt, startTime, endTime):
    '''
    DESCRIPTION
    Here we generate some Gaussian-like function discritized by the time-step 'dt'. More precisely, we create a gaussian
    like filter and then evaluate it at evenly spaced points and store this as a vector. This can then be used as one
    of the filters that we use to in some CIF.

    FUNCTIONAL FORM
    f(x) = a*exp(-(x-b)**2/(2*c)**2)

    :param a:
    :param b:
    :param c:
    :param startTime:
    :param stopTime:
    :param dt:
    :return:
    '''

    xvals = numpy.arange(start=startTime, stop=endTime, step=dt)
    convVector = a * numpy.exp(-1.*(((xvals-b)**2.)/(2.*c**2.)))
    return convVector


def halfGaussianFilterGenerator(a, b, c, dt, startTime, endTime):
    '''
    DESCRIPTION
    Here we generate some Gaussian-like function discritized by the time-step 'dt'. More precisely, we create a gaussian
    like filter and then evaluate it at evenly spaced points and store this as a vector. This can then be used as one
    of the filters that we use to in some CIF.

    The DIFFERNCE between this an your standard Gaussian, is that at input x>=0 we set the output to be = 0.

    '''

    xvals = numpy.arange(start=startTime, stop=endTime, step=dt)
    convVector = a * numpy.exp(-1. * (((xvals - b) ** 2.) / (2. * c ** 2.)))
    convVector[numpy.argwhere(xvals>=0)] = 0.
    return convVector




########################################################################################################################
# Convolutionl Functions
########################################################################################################################


def convolveSpikeTimes(convVector,spikeTimes,startTime,stopTime,dt):
    '''
    DESCRIPTION
    Given our spikeTime inputs, along with other bits, we calculate the values of the filter with the spike times.

    :param convVector:
    :param spikeTimes:
    :param startTime:
    :param stopTime:
    :param dt:
    :return:
    '''

    #timeArray = numpy.arange(start=startTime, stop=stopTime+dt, step=dt)
    boolArray = numpy.zeros(int((stopTime-startTime)/dt)+1) #We add one for the endpoint

    # Find indicies of the zeros array that we will want to set to a spike
    indsArray = ((spikeTimes - startTime)/dt)

    # Set the boolArray to have spikes where they should
    boolArray[indsArray.astype(int)] = 1.

    #Now convolve boolsArray with the convolution vector
    outputArray = numpy.convolve(boolArray, convVector, mode='full')

    #Output the entire array excpet the last values because they are not causel
    return outputArray[:-(convVector.size-1)]


########################################################################################################################
# Generation of Coefficients and Tie-Togethers
########################################################################################################################

# We need a set of algorithms to create and track what the coefficients we have are correlated too. MUCH of what will be
# here is perhaps overly hardcoded and needs to be broken up in future efforst.

def generateCifCompanionArrays(neuronGroups, neuronGroupStruc, neuronSynapses, filters):
    '''



    :param neuronGroups:
    :param neuronGroupStruc: is a list of numpy.ndarrays of the following form,
           numpy.array[[neuronId, x, y, z, groupId], ...]
    :param neuronSynapses: is a list of numpy.ndarrays of the following form,
           numpy.array[[neuronId-pre, neuronId-post, groupId-pre, groupId-post], ...]
    :param filters: is a list of lists of the following,
           where [neuronGroup, 'bias'] means we have a coefficient for each neuron in the group
           where [neuronGroup, 'gausLinSelf', a, b, c, timeLength, dt, timeDelays] means we add a linear filter for each neuron in said group tracking itself
           where [neuronGroup, 'gausLinSyn', neuronGroupTo, a, b, c, timeLength, dt, timeDelays] means we add linear summation of of synapses from
           where [neuronGroup, 'gausQuadSynSelf', neuronGroupTo, a, b, c, timeLength, dt, timeDelays] means we added quadratic coefficients from self comparison to input synapses
    :return:
    '''

    # Eumerate what we are building ####################################################################################
    # gausFiltersIdentities : numpy.ndarray, 2-D, each row has [a, b, c, timeLength, dt, time-delay, filterIndex]
    # gausFilters: numpy.ndarray, 2-D, each row has [a, b, c, timeLength, dt, truncation]
    # neuronFilters: numpy.ndarray, 2-D, each row has [neuron-group, neuron-index, filterIndex]


    # Build Filter Lists ##########################################################################################
    # We need to iterate through all the filter options and build a list or arrays of filters and a sub-array with the
    # actual filters. This is so we in the end perform the minimal amount of convolution operations.

    gausFiltersIdentities = numpy.array([[]])
    gausFilters = numpy.array([[]]) #

    # Iterate over the list of lists, each sub-list is a type of filter
    for i in range(len(filters)):

        if filters[i][1] == 'bias':
            newFilter = []
        elif filters[i][1] == 'gausLinSelf':
            newFilter = filters[i][2:]
        elif filters[i][1] == 'gausLinSyn':
            newFilter = filters[i][3:]
        elif filters[i][1] == 'gausQuadSynSelf':
            newFilter = filters[i][3:]
        else:
            print('\n\nERROR.... ERROR .... One of tye filter types cannot be identified ...\n\n')

        if len(newFilter)>0:
            # Parse over the different time-delay values to see if we need to add filters or truncated versions of said
            # filters.
            for ind, val in enumerate(newFilter[5]):
                if gausFiltersIdentities.size == 0:
                    # No filters have yet been added, so we need to make the first.
                    gausFiltersIdentities = numpy.array([newFilter[0], newFilter[1], newFilter[2], newFilter[3], newFilter[4], val, 0])
                    truncation = min(newFilter[3], 0.0-val)
                    gausFilters = numpy.array([newFilter[0], newFilter[1], newFilter[2], newFilter[3], newFilter[4], truncation])
                else:
                    # Check if our filter is already stored
                    indexF = numpy.where((gausFiltersIdentities[:,0:6]==(newFilter[0], newFilter[1], newFilter[2], newFilter[3], newFilter[4], val)).all(axis=1))
                    if indexF[0].size==0:
                        # Find out if the filter needs to be truncated
                        truncation = min(newFilter[3], 0.0 - val)
                        # Find out if the truncated filter is stored already
                        indexG = numpy.where((gausFilters==(newFilter[0], newFilter[1], newFilter[2], newFilter[3], newFilter[4], truncation)).all(axis=1))
                        if indexG[0].size==0:
                            #Both the filter and the truncated version need to be added!
                            gausFiltersIdentities = numpy.vstack((gausFiltersIdentities, numpy.array([newFilter[0], newFilter[1], newFilter[2], newFilter[3], newFilter[4], val, gausFilters.shape[0]])))
                            gausFilters = numpy.vstack((gausFilters, numpy.array([newFilter[0], newFilter[1], newFilter[2], newFilter[3], newFilter[4], truncation])))
                        else:
                            #Only the filter identity needs to be stored
                            gausFiltersIdentities = numpy.vstack((gausFiltersIdentities, numpy.array([newFilter[0], newFilter[1], newFilter[2], newFilter[3], newFilter[4], val, indexG[0][0]])))


    # Match Neurons to Filters #########################################################################################
    # Now we need to go through, and find which neurons need to be filtered, and denote by which filters

    neuronFilters = numpy.array([])

    # Iterate over the list of filters given as input
    for i in range(len(filters)):

        if filters[i][1] == 'gausLinSelf':
            # Iterate over the different time-delays for the given filter
            for ind, val in enumerate(filters[i][7]):
                # Find the appropriate filterIndex
                indexF = numpy.where((gausFiltersIdentities[:,0:6]==(filters[i][2], filters[i][3], filters[i][4], filters[i][5], filters[i][6], val)).all(axis=1))
                filterIndex = gausFiltersIdentities[indexF[0][0],6].copy() # Do I want to use the copy() to protect?
                # Iterate over the group of neurons
                for j in range(neuronGroupStruc[filters[i][0]].shape[0]):
                    indexG = numpy.where((neuronFilters==(filters[i][0], j, filterIndex)).all(axis=1))
                    if indexG[0].size==0:
                        neuronFilters = numpy.vstack((neuronFilters,numpy.array([filters[i][0], j, filterIndex])))
        elif filters[i][1] == 'gausLinSyn':
            # Iterate overr the different time-dealys for the given filter
            for ind, val in enumerate(filters[i][8]):
                # Find the appropriate filterIndex
                indexF = numpy.where((gausFiltersIdentities[:, 0:6] == (filters[i][3], filters[i][4], filters[i][5], filters[i][6], filters[i][7], val)).all(axis=1))
                filterIndex = gausFiltersIdentities[indexF[0][0], 6].copy()  # Do I want to use the copy() to protect?
                # Now I want to look at the neurons I care about! In order to do this I must find the synapses with
                # the correct groups.
                for synArray in neuronSynapses:
                    if synArray[0,2]==filters[i][2] and  synArray[0,3]==filters[i][0]:
                        for unNeu in numpy.unique(synArray[:,0]):
                            indexG = numpy.where((neuronFilters==(filters[i][2], unNeu, filterIndex)).all(axis=1))
                            if indexG[0].size==0:
                                neuronFilters = numpy.vstack((neuronFilters,numpy.array([filters[i,2], unNeu, filterIndex])))
        elif filters[i][1] == 'gausQuadSynSelf':
            # Iterate over the different time-delays for the given filter
            for ind, val in enumerate(filters[i][8]):
                # Find the appropriate filterIndex
                indexF = numpy.where((gausFiltersIdentities[:, 0:6] == (filters[i][3], filters[i][4], filters[i][5], filters[i][6], filters[i][7], val)).all(axis=1))
                filterIndex = gausFiltersIdentities[indexF[0][0], 6].copy()  # Do I want to use the copy() to protect?
                # First I want to iterate over the group of my neurons
                for j in range(neuronGroupStruc[filters[i][0]].shape[0]):
                    indexG = numpy.where((neuronFilters==(filters[i][0], j, filterIndex)).all(axis=1))
                    if indexG[0].size==0:
                        neuronFilters = numpy.vstack((neuronFilters,numpy.array([filters[i][0], j, filterIndex])))
                # Secondly I want to look at the neurons that connect into my neurons
                for synArray in neuronSynapses:
                    if synArray[0,2]==filters[i][2] and synArray[0,3]==filters[i][0]:
                        for unNeu in numpy.unique(synArray[:, 0]):
                            indexG = numpy.where((neuronFilters == (filters[i][2], unNeu, filterIndex)).all(axis=1))
                            if indexG[0].size == 0:
                                neuronFilters = numpy.vstack((neuronFilters, numpy.array([filters[i, 2], unNeu, filterIndex])))


                # Now I want to look at the neurons I care about! In order to do this I must find the synapses with
                # the correct groups.






    # Make Filter Convolution Kernels ##################################################################################
    # Now we make the actual kernels that we convolve the spike time series with


    # Create Filter-Neuron Array Listing ###############################################################################
    # This will tell the user or the code where to find a specific neurons


















########################################################################################################################
# Internal Unit Testing
########################################################################################################################

if __name__ == "__main__":

    #Test function generation
    vectorVals = gaussianFilterGenerator(a=1.0, b=0., c=.01, dt=.001, startTime=-.1, endTime=.1)
    plt.figure(1)
    plt.plot(vectorVals)
    #vectorVals = halfGaussianFilterGenerator(a=1.0, b=0., c=.01, dt=.001, startTime=-.1, endTime=.1)
    #plt.plot(vectorVals)

    #Test the convolution to spike times thing
    spikeTimes = numpy.array([1.001, 1.012, 1.045, 1.09, 1.231, 1.298, 1.521, 1.529, 1.545])
    evaled = convolveSpikeTimes(vectorVals, spikeTimes, startTime=1.000, stopTime=1.6, dt=.001)
    startTime = 1.6 - (evaled.size-1)*.001
    print(startTime)
    timesArray = numpy.arange(start=startTime, stop=1.6+.001, step=.001)
    print(timesArray.shape)
    plt.figure(2)
    plt.plot(timesArray, evaled)
    plt.show()
