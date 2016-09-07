# A.Lonsberry
# June 2016
#
# DESCRIPTION
#   Here all the different elements are put together and run
import time
import math
import numpy
import matplotlib.pyplot as plt
from brian2 import *

import condIntenFunc_0_0 #Input-spikes generations function
import condIntenFunc_0_1 #Network organization functions
import condIntenFunc_0_2 #Network simulation functions
import condIntenFunc_0_3 #CIF Functions

if __name__ == "__main__":

    ####################################################################################################################
    # DEFINE NETWORK
    ####################################################################################################################

    # Here we want one neuron that stimulates a larger network with some predefiend stimulus
    #           {    O - O - O   O - O -- O   }
    #           {   /    |    \ /    |   /|   }
    #   O-------{--O --- O --- O     |  / O   }
    #    \      {       /|    / \    | /      }
    #     ------{----O - O ----- O - O ---O   }
    #  INPUT             MAIN-NETWORK

    # Input layer/body-of-neurons: a single neuron
    inputNetStruc = condIntenFunc_0_1.createNetLayer(xN=2, groupId=0)

    # Main layer/body-of-neurons: a large grid of neurons
    mainNetStruc = condIntenFunc_0_1.createNetLayer(xN=8, yN=8, groupId=1)

    # The connection structure from input-neurons to main-neurons is defined by probability to each layer. We here have
    # to define this function.
    def layerConnProbFunc(x1, y1, z1, x2, y2, z2):
        output = .1
        return output

    # Define connection probabilities from input-body to main-body
    inputToMainConnProb = condIntenFunc_0_1.createConnectionProb(inputNetStruc, mainNetStruc, layerConnProbFunc)

    # Define connections from input-body to main-body
    inputToMainConnections = condIntenFunc_0_1.createConnectionArray(inputToMainConnProb, groupId1=0, groupId2=1)

    # The connection structure between the mian-boud of neuorns is defiend by a probability function which we have to
    # here define
    def bodyConnProbFunc(x1, y1, z1, x2, y2, z2):
        output = 0.30
        return output

    # Def ine connection probabilites from main-body neurons to main-body neurons
    mainToMainConnProb = condIntenFunc_0_1.createConnectionProb(mainNetStruc, mainNetStruc, bodyConnProbFunc )

    # Define connection from main-body to main-body
    mainToMainConnections = condIntenFunc_0_1.createConnectionArray(mainToMainConnProb, groupId1=1, groupId2=1)

    ####################################################################################################################
    # CREATE INPUT DATA
    ####################################################################################################################

    # Assuming 1 single input neuron, we are going to create data that it inputs into our main network. Here we are
    # going to create an input of alternating firing rates.

    inputArray = condIntenFunc_0_0.multiPoissonStreams(2, startTime=0.0, stopTime=100.0, rateBounds=[20., 40.], intervalBounds=[2., 8.])
    print(inputArray)

    ####################################################################################################################
    # RUN SIMULATION WITH THE ADAPTING SYNAPTIC WEIGHTS
    ####################################################################################################################

    simResultsList = condIntenFunc_0_2.singleLayerSimulationNewAdapt(inputNetStruc, mainNetStruc, inputToMainConnections,
                                                                     numpy.array([]), mainToMainConnections,
                                                                     numpy.array([]), inputArray)

    ####################################################################################################################
    # CHECK RESULTS (Debugging)
    ####################################################################################################################

    # This is primarily for debugging. I want to make sure it looks like the simulations are working as expected.
    M0 = simResultsList[0] # input-layer neuron(s) spike monitor
    M1 = simResultsList[1] # main-body neuron(s) spike monitor
    M2 = simResultsList[2] # main-body neuron(s) spike-rate monitor
    S  = simResultsList[3] #
    S2 = simResultsList[4]
    M3 = simResultsList[5]
    M4 = simResultsList[6]
    M5 = simResultsList[7]

    """
    # Plot main-body neuron group smooth firing rate.
    plt.figure(1)
    plt.plot(M2.t, M2.smooth_rate(window='gaussian', width=3. * second) / Hz, '-b')
    plt.xlabel('Time')
    plt.ylabel('Firing Rate (Hz)')
    plt.title('Main Neuron Group Avg. Firing Rate vs. Time')

    # Plot Raster plot of
    plt.figure(2)
    plt.subplot(211)
    plt.plot(M0.t, M0.i, '.k')
    plt.xlabel('Time')
    plt.subplot(212)
    plt.plot(M1.t, M1.i, '.k')
    plt.xlabel('Time')

    # Plot noisey input and voltage
    plt.figure(3)
    plt.subplot(211)
    plt.plot(M3.t, M3.I[0])
    plt.subplot(212)
    plt.plot(M3.t, M3.v[0])

    # Plot Synaptic weight changes
    plt.figure(4)
    plt.plot(M4.t, M4.w[0])

    #plt.show()
    """

    ####################################################################################################################
    # RUN NETWORK WITHOUT ADAPTING SYNAPSES
    ####################################################################################################################
    # This is done such that data can be collected in order to calculate and find conditional intensity functions.

    # Find resultant synaptic weights. Firstly make an empty matrix and then fill with the appropriate weights.
    inputToMainConnWeights = numpy.zeros(inputToMainConnections.shape[0])
    for i in range(inputToMainConnections.shape[0]):
        inputToMainConnWeights[i] = S.w[i]
    mainToMainConnWeights = numpy.zeros(mainToMainConnections.shape[0])
    for i in range(mainToMainConnections.shape[0]):
        mainToMainConnWeights[i] = S2.w[i]

    # Now we need to rerun with network that does
    simResultsList = condIntenFunc_0_2.singleLayerSimulationNew(inputNetStruc, mainNetStruc,
                                                                inputToMainConnections,
                                                                inputToMainConnWeights, mainToMainConnections,
                                                                mainToMainConnWeights, inputArray)

    # This is primarily for debugging. I want to make sure it looks like the simulations are working as expected.
    M00 = simResultsList[0]  # input-layer neuron(s) spike monitor
    M11 = simResultsList[1]  # main-body neuron(s) spike monitor
    M22 = simResultsList[2]  # main-body neuron(s) spike-rate monitor
    S = simResultsList[3]  #
    S2 = simResultsList[4]
    M33 = simResultsList[5]
    M44 = simResultsList[6]
    M55 = simResultsList[7]


    # Plot main-body neuron group smooth firing rate.
    plt.figure(5)
    plt.plot(M22.t, M22.smooth_rate(window='gaussian', width=3. * second) / Hz, '-b')
    plt.xlabel('Time')
    plt.ylabel('Firing Rate (Hz)')
    plt.title('Main Neuron Group Avg. Firing Rate vs. Time')

    # Plot Raster plot of
    plt.figure(6)
    plt.subplot(211)
    plt.plot(M00.t, M00.i, '.k')
    plt.xlabel('Time')
    plt.subplot(212)
    plt.plot(M11.t, M11.i, '.k')
    plt.xlabel('Time')

    # Plot noisey input and voltage
    plt.figure(7)
    plt.subplot(211)
    plt.plot(M33.t, M33.I[0])
    plt.subplot(212)
    plt.plot(M33.t, M33.v[0])

    # Plot Synaptic weight changes
    plt.figure(8)
    plt.plot(M44.t, M44.w[0])

    plt.show()

    ####################################################################################################################
    # CREATE CIF FUNCTIONS
    ####################################################################################################################

    # Create a new conditional-intensity-function manager object
    myFilterMan = condIntenFunc_0_3.condItenFuncManager(.001)  # We have to give it the dt we are using
    # We only care about doing this for a couple of neurons at most at the moment
    # neuronGroup0Ids = numpy.arange(0, inputNetStruc.shape[0], 1)
    neuronGroup1Ids = numpy.arange(0, 4, 1)
    gausFiltDelays = numpy.array([0., .005, .01, .015, .02, .025, .03, .035, .040, .045, .05, .055, .065, .075, .08, .085, .09, .095, 0.1, .105, .11, .115, .12, .125, .13, .135, .140, .145, .15, .155, .165, .175, .18, .185, .19, .195])
    inputFunctionType = ['bias', 1, neuronGroup1Ids]
    myFilterMan.addNewFilter(inputFunctionType)
    inputFunctionType = ['gausLinSelf', 1, neuronGroup1Ids, 1., .0, .002, .1, .001, gausFiltDelays]
    myFilterMan.addNewFilter(inputFunctionType)
    inputFunctionType = ['gausLinSyn', 1, neuronGroup1Ids, mainToMainConnections, 1., .0, .002, .1, .001, gausFiltDelays]
    myFilterMan.addNewFilter(inputFunctionType)
    inputFunctionType = ['gausLinSyn', 1, neuronGroup1Ids, inputToMainConnections, 1., .0, .002, .1, .001, gausFiltDelays]
    myFilterMan.addNewFilter(inputFunctionType)
    myFilterMan.filterManager.initFilterRecorders()
    # Now actually give the data to be presented, note we need to get ride of the units brian2 has in the spike monitors
    # we do this by creating the following objects below.
    class someSpikeObj:
        def __init__(self, inds, times):
            self.t = times
            self.i = inds
        def correctTimes(self, startTime, endTime):
            timeCorrectedInds = numpy.where( (self.t>=startTime) & (self.t<=endTime))
            self.t = self.t[timeCorrectedInds]
            self.i = self.i[timeCorrectedInds]
            print(type(self.t))
            print(self.t)
    # Now we want to try and train using gradient descent! We wrap this into a fucntion here below
    def testSgdAndPrections(myFilterMan, trainIterations, cifObjNumb):
        # Grab the cif-object
        cifObj = myFilterMan.condIntenFunctions[cifObjNumb]
        # Iteritively train using SGD method
        for i in range(trainIterations):
            cifObj.updateCifCoeffsSGD(.01)

    # We have a lot of data, we want to break it up into smaller pieces and run SGD on the data.
    #timeBounds = numpy.arange(0, 100, 5)
    timeBounds = numpy.arange(0, 100, 5)
    for iterOver in range(10):
        for i in range(timeBounds.shape[0]-1):
            timeStr = "Training from time = " + str(timeBounds[i]) + " to time = " + str(timeBounds[i+1])
            print(timeStr)
            spikObj0 = someSpikeObj(numpy.asarray(M00.i), numpy.asarray(M00.t))
            spikObj0.correctTimes(timeBounds[i], timeBounds[i+1])
            spikObj1 = someSpikeObj(numpy.asarray(M11.i), numpy.asarray(M11.t))
            spikObj1.correctTimes(timeBounds[i], timeBounds[i + 1])
            t0 = time.clock()
            myFilterMan.parseNewSpikeTrainData([spikObj0, spikObj1], [0, 1], startTime=timeBounds[i], endTime=timeBounds[i+1])
            rTime = time.clock() - t0
            rStr = "...\nTime to parse new data = " + str(rTime)
            print(rStr)
            t0 = time.clock()
            myFilterMan.prepareParsedData()
            rTime = time.clock() - t0
            rStr = "Time to prepare data = " + str(rTime)
            print(rStr)
            # Run the test method we just created on one particular neuron in question
            t0 = time.clock()
            testSgdAndPrections(myFilterMan, 2, 0)
            rTime = time.clock() - t0
            rStr = "Time to SGD data = " + str(rTime)
            print(rStr)
            # Clear out residual data
            myFilterMan.garbageCollect(endTime=timeBounds[i+1])
        # Clear out all memory and reset
        myFilterMan.clearMemory()



    # Now we want to see how well the SGD worked
    plt.figure(1)
    spikObj0 = someSpikeObj(numpy.asarray(M00.i), numpy.asarray(M00.t))
    spikObj0.correctTimes(0, 40)
    spikObj1 = someSpikeObj(numpy.asarray(M11.i), numpy.asarray(M11.t))
    spikObj1.correctTimes(0, 40)
    myFilterMan.parseNewSpikeTrainData([spikObj0, spikObj1], [0, 1], startTime=0, endTime=40)
    myFilterMan.prepareParsedData()
    cifObj = myFilterMan.condIntenFunctions[0]
    plt.plot(cifObj.spikeTimeArray[-cifObj.dataMatrix.shape[0]:], cifObj.spikeArray[-cifObj.dataMatrix.shape[0]:], '--', linewidth=3.0)
    # Generate initial predictions, and plot them
    probs = cifObj.generateSpikeProbabilities()
    plt.plot(cifObj.spikeTimeArray[-probs.size:], probs, linewidth=4.0)
    print(cifObj.spikeArray)
    print(cifObj.spikeTimeArray)
    print(M11.i.shape)
    print(M11.t.shape)
    theInds = numpy.where(M11.i==0)
    print(type(M11.t))
    print(M11.t[theInds[0]])
    plt.show()



