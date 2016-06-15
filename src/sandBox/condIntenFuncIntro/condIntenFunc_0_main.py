# A.Lonsberry
# June 2016
#
# DESCRIPTION
#   Here all the different elements are put together and run

import math
import numpy
import matplotlib.pyplot as plt
from brian2 import *

import condIntenFunc_0_0 #Input-spikes generations function
import condIntenFunc_0_1 #Network organization functions
import condIntenFunc_0_2 #Network simulation functions

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
    mainNetStruc = condIntenFunc_0_1.createNetLayer(xN=10, yN=10, groupId=1)

    # The connection structure from input-neurons to main-neurons is defined by probability to each layer. We here have
    # to define this function.
    def layerConnProbFunc(x1, y1, z1, x2, y2, z2):
        output = .03
        return output

    # Define connection probabilities from input-body to main-body
    inputToMainConnProb = condIntenFunc_0_1.createConnectionProb(inputNetStruc, mainNetStruc, layerConnProbFunc)

    # Define connections from input-body to main-body
    inputToMainConnections = condIntenFunc_0_1.createConnectionArray(inputToMainConnProb, groupId1=0, groupId2=1)

    # The connection structure between the mian-boud of neuorns is defiend by a probability function which we have to
    # here define
    def bodyConnProbFunc(x1, y1, z1, x2, y2, z2):
        output = 0.20
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

    inputArray = condIntenFunc_0_0.multiPoissonStreams(2, startTime=0.0, stopTime=5000.0, rateBounds=[10., 20.], intervalBounds=[2., 8.])
    print(inputArray)

    ####################################################################################################################
    # Run Simulation
    ####################################################################################################################

    simResultsList = condIntenFunc_0_2.singleLayerSimulationNew(inputNetStruc, mainNetStruc, inputToMainConnections,
                                                                numpy.array([]), mainToMainConnections,
                                                                numpy.array([]), inputArray)


    ####################################################################################################################
    # Check results
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

    plt.show()

