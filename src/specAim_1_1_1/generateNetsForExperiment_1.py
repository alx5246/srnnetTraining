# A.Lonsberry
# November 2016
#
# DESCRIPTION
#   Here I am going to have methods to create, test, and alter (evolutionary/genetic) networks that will be used in the
#   experiments. To be very clear, here I will essentially create and save the networks. These networks can then be
#   loaded into a script to do all the experiments.
#

import matplotlib.pyplot as plt
import numpy

# I had a difficult time getting things to import, so now I have however, go about fixing this by going up a directoy
# which I think I can do in the following way.
import os
import sys
import pickle
sys.path.append(os.path.abspath("../"))

# Import some of my other moduels in my repository
from networkGenerators import networkGenerator


def genProtocal_1_NetStrcut(numInputs, numOutputs, saveName=[]):
    """
    DESCRIPTION
    This will create the appropriate networks for first experiemts/protocol_1. Our networks here will feature multiple
    inputs and 1 output neuron (most generally). This method is pretty simple as it generates all the network structure
    but does not add any network dynamics or models. This method returns and can optionally also save data as well.
    :param numInputs: integer, the number of input signals
    :param numOutputs: integer, the number of mian-layer signals (should be 1)
    :return: [inputNetStruc, mainNetStruct, inputToMainConnections]
    """
    # Create the input layer, which is simply the input signals in spike train form
    inputNetStruc = networkGenerator.createNetLayer(xN=numInputs, groupId=0)
    # Create the output layer,
    mainNetStruct = networkGenerator.createNetLayer(xN=numOutputs, groupId=1)
    # Define the method that controls the input-layer to final-layer connection prob
    def layerConProbFunc(x, y, z, x1, y1, z1):
        output = 1.
        return output
    # Define connection probabilities from input to main-layer
    inputToMainConnProb = networkGenerator.createConnectionProb(inputNetStruc,mainNetStruct,layerConProbFunc)
    # Define connections from input to main
    inputToMainConnections = networkGenerator.createConnectionArray(inputToMainConnProb, groupId1=0, groupId2=1)

    # Now we save the network to a file (only if a file-name is given)
    if isinstance(saveName,str):
        if saveName[-4:]=='.pkl':
            # Use pickle to save the data we generated above
            outputList = [inputNetStruc, mainNetStruct, inputToMainConnections]  # This is what we will be saving
            outputFile = open(saveName, "wb")
            pickle.dump(outputList, outputFile)
            outputFile.close()
            # We also want to generate a description file automatically as well
            descpFName = saveName[:-4] + 'DESC.txt' # The file name of the description files
            file = open(descpFName, 'w')
            file.write("Experiment/Protocol 1 Network Structure for file...")
            file.write("\n..."+saveName+"\n")
            file.write("\nData here is generated from the genProtocal_1_NetStrut\n")
            file.write("\nNumber of inputs........."+str(numInputs))
            file.write("\nNumber of outputs........"+str(numOutputs))
            file.write("\nConnectivity = 100%")

    return [inputNetStruc, mainNetStruct, inputToMainConnections]


def singleLayerSim_00(inputNetStruc, mainNetStruc, inputToMainConnections, inputToMainConnWeights, mainToMainConnections, mainToMainConnWeights, inputSpikes):
    '''
    DESCRIPTION
    This runs a new simulation with a Izh. RS neurons (slightly modified to include regular current). I do not have any
    STDP or Homeostatsis developed in Carlson et al. 2013. Here I will
    :param inputNetStruc: structure of input layer of neurons 2D numpy.array([[]]), size = numb. neurons X 4, in each
           row we have [the neuron index, x-position, y-position, z-position]
    :param mainNetStruc: structure of input layer of neurons 2D numpy.array([[]]), size = numb. neurons X 4, in each
           row we have [the neuron index, x-position, y-position, z-position]
    :param inputToMainConnections: numpy.ndarray, nX4, in each row [[presynaptic Neuron ID, postsynaptic neuron ID, presynatpic neuron group,
           postsynaptic neuron group] ...]
    :param mainToMainConnections: numpy.ndarray, nX4, in each row [[presynaptic Neuron ID, postsynaptic neuron ID, presynatpic neuron group,
           postsynaptic neuron group] ...]
    :param inputSpikes: 2D numpy.array, number of rows equal total number of spikes, there are two columns. In the first column is
           the index of the neuron that spiked, and in the second column is the time of said spike. The array is sorted by
           time in asscending order. This input will be given to our input layer.
    :return:
    '''

########################################################################################################################
# RUNNING METHODS
########################################################################################################################
# Some the methods will be used to create acutal data and some will be used for testing. The differentiation will be
# given.

if __name__ == "__main__":

    print('\n\nRUNNING ')

    ####################################################################################################################
    # TESTING METHODS
    ####################################################################################################################
    # You have to uncomment the method calls to run!
    genProtocal_1_NetStrcut(numInputs=10, numOutputs=1, saveName='networks/testNet.pkl')



