# A.Lonsberry
# November 2016
#
# DESCRIPTION
#   Here I am going to implement methods to alter the format of data.

import matplotlib.pyplot as plt
import numpy

# I had a difficult time getting things to import, so now I have however, go about fixing this by going up a directoy
# which I think I can do in the following way.
import os
import sys
import pickle
sys.path.append(os.path.abspath("../"))


def convertSpikeListsToBrainInput(spikeLists):
    """
    DESCRIPTION
    The methods in spikeDataGenerators/analogToSpikes.py have methods that produce spikes in list form. The typical
    spike output is a list of lists, where each sublist has an ordered set of spike times from a neuron or input signal.
    Here we want to convert this into a form typically used by brian module, where input in a nX2 numpy.array with the
    first column having the neuron id, and the second column having teh spike time. This list needs to be sorted by
    time in asscending order.
    :param spikeLists: list of list [ [.2, .4, 19.1], .... ] where each sublist is the spike times for a neuron or input
           signal
    :return: nX2 numpy.array(), where the first column has the signal ID, the second column has spike time, and the
             whole thing is  sorted by time in asscending order.
    """

    # First I have function that will take in a list of lists and convert to correct numpy.array() structure
    def altSpikeStruct(listOfSpikeLists):
        # Find how many spikes in total we will have
        numbInds = 0
        for someList in listOfSpikeLists:
            numbInds += len(someList)
        # Create the full numpy array, and then fill it
        outArray = numpy.zeros((numbInds, 2))
        index = 0
        for i,someList in enumerate(listOfSpikeLists):
            outArray[index:index+len(someList), 0] = i
            outArray[index:index+len(someList), 1] = numpy.asarray(someList)
            index += len(someList)
        return outArray

    # Now I will define a function so I can recurse into the spikeLists thing
    def incrementInside(inputList):
        if isinstance(inputList[0],list):
            if isinstance(inputList[0][0],list):
                # If we have a list within a list, we are not yet deep enough, though we will have to begin to iterate
                subList = []
                for i in range(len(inputList)):
                    subList.append(incrementInside(inputList[0]))
                return subList
            else:
                # We are deep enough
                return altSpikeStruct(inputList)
        else:
            # If we are deep enough we call the spike-times to numpy.array function
            return altSpikeStruct([inputList])

    # Now we will apply the recursion
    output = incrementInside(spikeLists)
    return output




########################################################################################################################
# TESTING METHDOS
########################################################################################################################
# Some the methods win.

if __name__ == "__main__":

    import random

    # TESTING convertSpikeListsToBrainInput()
    # Test 1 : a spike list of spike times
    spikeList = []
    for i in range(5):
        spikeList.append(random.random())
        if len(spikeList)>2:
            spikeList[-1] = spikeList[-1] + spikeList[-2]
    outputSpikesFormatted = convertSpikeListsToBrainInput(spikeList)
    print("\n")
    print(outputSpikesFormatted )

    # Test 2: two spike lists in one list
    spikeList = []
    for i in range(2):
        subSpikeList = []
        for j in range(5):
            subSpikeList.append(random.random())
            if len(subSpikeList)>2:
                subSpikeList[-1] = subSpikeList[-1] + subSpikeList[-2]
        spikeList.append(subSpikeList)
    outputSpikesFormatted = convertSpikeListsToBrainInput(spikeList)
    print("\n")
    print(outputSpikesFormatted)

    # Test 3: Two spike lists in a bigger list
    spikeList = [[spikeList, spikeList],[ spikeList]]
    outputSpikesFormatted = convertSpikeListsToBrainInput(spikeList)
    print("\n")
    print(outputSpikesFormatted)

