# A.Lonsberry
# June 2016
#
# DESCRIPTION
#   Here are functions to create and annotate some network, this is done to keep track of a network and its configuration
#   if we need to reinstantiate it later.

import numpy
import random

def createNetLayer(xN, yN=0, zN=0, groupId=0):
    '''
    DESCRIPTION
    Creates an array of neuron positions inside a cube lattice structure. This tells us the network position of a
    particular neuron.
    :param xN: integer, the number of neurons along the x-axis of some lattice of neurons
    :param yN: integer, the number of neurons along the y-axis of some lattice of neurons
    :param zN: integer, the number of neurons along the z-axis of some lattice of neurons
    :param groupId: int, tells us which neuron group we are claiming
    :return: 2D numpy.ndarray, size = numb. neurons X 4, in each row we have [the neuron index, x-position, y-position, z-position]
    '''
    totalNumbNeurons = max(xN*yN*zN, xN*yN, yN*zN, zN*xN, xN, yN, zN)
    # Perform some error handling and number alterations
    xN = int(xN)
    yN = int(yN)
    zN = int(zN)
    groupId = int(groupId)
    if xN <= 0:
        xN = 1
    if yN <= 0:
        yN = 1
    if zN <= 0:
        zN = 1
    # Creat the output, which is a numpy.array which lists where each neuron is within a 3D structure, and what its ID
    # is in that structure (the neuorn index). Then fill in this array.
    netNeurons = numpy.zeros((totalNumbNeurons, 5), dtype=int)
    nCount = 0
    for i in range(xN):
        for j in range(yN):
            for k in range(zN):
                netNeurons[nCount, 0] = nCount
                netNeurons[nCount, 1] = i
                netNeurons[nCount, 2] = j
                netNeurons[nCount, 3] = k
                netNeurons[nCount, 4] = groupId
                nCount += 1
    return netNeurons


def createConnectionProb(net1, net2, probGenerator):
    '''
    DESCRIPTION
    Given some function (probGenerator) and two networks (net1 and net2), this thing will produce an array of the
    probabilities of connections between neurons in each network.
    :param net1: 2D numpy.array, with [ [neuronID, x-position, y-position, z-postion]...]
    :param net2: 2D numpy.array, with [ [neuronID, x-position, y-position, z-postion]...]
    :param probGenerator: a python function that takes as input [x, y, z, x1, y1, z1] and returns probability of
           connection from net1 neurons (x, y, z) to net2 neurons (x1, y1, z1).
    :return: 2D numpy.ndarray of size net1.shape[0] X net2.shape[0], values over [0., 1.] (given input function)
    '''
    # Initialize the output array, which will hold the probabilities of connection from neuron i (row) to neuron j (column)
    connectionProb = numpy.zeros((net1.shape[0], net2.shape[0]))
    for i in range(net1.shape[0]):
        for j in range(net2.shape[0]):
            connectionProb[i, j] = probGenerator(net1[i,1], net1[i,2], net1[i,3], net2[j,1], net2[j,2], net2[j,3])
    return connectionProb


def createConnectionArray(layerProb, groupId1, groupId2):
    '''
    DESCRIPTION
    Given some connection probabilities, we will actually generate the connections and output what those connections
    actually are.
    :param layerProb: 2D numpy.array, nXm where n = number of neurons in groupId1, and m = number of neurons in groupId2
        where each index stores the probability of connection from neuron i (row) to neuron j (column)
    :param groupId1: int, the group the presynaptic neuron is within
    :param groupId2: int, the group the postsynatic neuron is within
    :return: numpy.ndarray, nX4, in each row [[presynaptic Neuron ID, postsynaptic neuron ID, presynatpic neuron group,
        postsynaptic neuron group] ...]. This numpy.array is how we will store connections between neuron groups.
    '''
    # Inital array where we store the connections bits.
    connectionArray = numpy.zeros((layerProb.shape[0], layerProb.shape[1]), dtype=int)
    # Fill connectionArray with ones on a probabilistic basis
    for i in range(layerProb.shape[0]):
        for j in range(layerProb.shape[1]):
            if random.random() <= layerProb[i, j]:
                connectionArray[i, j] = 1
    #Now we want to trim down into an array with only pre- and post-synaptic connections
    trimmedArray = numpy.zeros((numpy.sum(connectionArray), 4), dtype=int)
    counter = 0
    for i in range(connectionArray.shape[0]):
        for j in range(connectionArray.shape[1]):
            if connectionArray[i, j]==1:
                trimmedArray[counter, 0] = i
                trimmedArray[counter, 1] = j
                trimmedArray[counter, 2] = groupId1
                trimmedArray[counter, 3] = groupId2
                counter += 1
    return trimmedArray


########################################################################################################################
# Internal Testing
########################################################################################################################
# I want to do some very simple unit testing on the functions above.

if __name__ == "__main__":

    print("Entering Unit Test ... ... \n")

    print("Making network ... ")
    print("   Printing out network ...")
    myFirstNet = createNetLayer(2, 0, 0)
    print(    myFirstNet)
    print("\n")

    print("Making second network ...")
    print("   Printing out sescond network")
    mySecondNet = createNetLayer(2, 2, 0)
    print(  mySecondNet)
    print("\n")

    print("Defining my first connection function ... ")
    def layerProb(x1, y1, z1, x2, y2, z2):
        if x2==0:
            output = 1.0
        else:
            output = 0.0
        return output
    print("Creating connection functions ...")
    net1ToNet2ConnProb = createConnectionProb(myFirstNet, mySecondNet, layerProb)
    print("   Printing connection prob ...")
    print(net1ToNet2ConnProb)
    print("\n")

    print("Creating connection array ...")
    net1ToNet2Conn = createConnectionArray(net1ToNet2ConnProb)
    print("   Pringint connection array ...")
    print(net1ToNet2Conn)
    print("\n")

    print("Defining my second connection function ...")
    def layerProb1(x1, y1, z1, x2, y2, z2):
        if x1 != x2 and y1 != y2 and z1 != z2:
            output = .5
        else:
            output = 0.
        return output
    net2ToNet2ConnProb = createConnectionProb(mySecondNet, mySecondNet, layerProb1)





