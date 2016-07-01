# S. Pickard
# June 2016
# (Using python 3.4.4, and brian2 2.0rc)
#
# DESCRIPTION
#   Function/method that produces an output
#   This analysis was inspired by Lazar et al.(2009) Figure 3.H
#   Goal is

import numpy as np

def spikeMon_To_Matrix(spikeTimeArray, NeurIndexArray):
    """
    FUNCTION DESCRIPTION
        This function takes as input two 1-D numpy.arrays of spike times (spikemon.t) and the corresponding neuron
        index array for each neuron index (spikemon.i). This function will generate a N x T array where each row
        corresponds to a neuron in the network (N dimensions) by the spike times (T dimensions).

    :param spikeTimeArray: n-D numpy.array, units are seconds, neuron spike times for each neuron stored in a numpy.array
    :param NeurIndexArray: n-D numpy.array, units are scalar, neuron index for each neuron stored in a numpy.array

    """

    #Toy Data (comment out when importing simulation data)
    # NeurIndexArray = [1,1,3,1,2,3,0]
    # spikeTimeArray = [.5,.8,1.1,1.6,2.5,2.9, 2.9]
    # print(spikeTimeArray)

    #Make input arrays into numpy arrays
    spikeTimeArray = np.array(spikeTimeArray)
    # print(spikeTimeArray)
    NeurIndexArray = np.array(NeurIndexArray)
    # print(NeurIndexArray)


    #Set parameters for loop
    iter = np.amax(NeurIndexArray)
    # print (iter)

    i = 0           #Counter for loop iteration tracking
    ind = 0         #tracker for neuron index positions
    NeurInd = []    #Holds positions for each neuron

    #For loop to make N number of arrays that hold the element positions corresponding to neuron firing time
    for i in range(iter+1):
        ind = np.where(NeurIndexArray == i)     #holder array; i value will correspond to neuron index; find all elements for each value of i
        NeurInd.append(ind)     #array holding element positions for each neuron; in sequential order


    # print ('NeurInd: ', NeurInd)
    # print ('NeurInd length: ', len(NeurInd))
    # print('NeurInd[0]: ', NeurInd[0])

    # print (NeurInd[0])
    # print (NeurInd[1])

    #Parameters for second loop
    j = 0               #Counter for loop iteration tracking
    NFT = 0             #temporary holder for fire times
    NeurFireTime = []   #Array of fire times in order of neuron index

    #For loop used to extract firing times corresponding to neuron index
    for j in range(iter+1):
        NFT = spikeTimeArray[NeurInd[j]]    #value holder; holds spike times corresponding to elements to corresponding neurons
        NeurFireTime.append(NFT)


    # print ('NeurFireTime: ', NeurFireTime)
    # # print('NeurFireTime length: ', len(NeurFireTime))
    # # print ('NeurFireTime[0]: ', NeurFireTime[0])
    #
    # # print (type(NeurFireTime))
    # # print (NeurFireTime.shape)
    #
    return (NeurFireTime)