# A.Lonsberry
# June 2016
#
# DESCRIPTION
#   Here all the different elements are put together and run

import math

import condIntenFunc_0_0 #Input-spikes generations function
import condIntenFunc_0_1 #Network organization functions

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
    #  INPUT             MAIN NETWORK

    # Input layer/body-of-neurons: a single neuron
    inputNetStruc = condIntenFunc_0_1.createNetLayer(xN=1)

    # Main layer/body-of-neurons: a large grid of neurons
    mainNetStruc = condIntenFunc_0_1.createNetLayer(xN=10, yN=10)

    # The connection structure from input-neurons to main-neurons is defined by probability to each layer. We here have
    # to define this function.
    def layerConnProbFunc(x1, y1, z1, x2, y2, z2):
        output = math.exp(-(x2))
        return output

    # Define connection probabilities from input-body to main-body
    inputToMainConnProb = condIntenFunc_0_1.createConnectionProb(inputNetStruc, mainNetStruc, layerConnProbFunc)
    # Define connections from input-body to main-body
    inputToMainConnections = condIntenFunc_0_1.createConnectionArray(inputToMainConnProb)


    ####################################################################################################################
    # CREATE INPUT DATA
    ####################################################################################################################

    #Assuming 1 single input neuron, we are going to create data that it inputs into our main network. Here we are going
    #to create an input of alternating firing rates.
    #spikeArray = condIntenFunc_0_0.poissonSpikeGen([5., 10., 5., 10., 5., 10, 5., 10.], [0., 10., 20., 30., 40., 50., 60., 70., 80.])

    #Now we have to alter the spikeArray so that it works with