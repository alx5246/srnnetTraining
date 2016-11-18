# A.Lonsberry
# November 2016
#
# DESCRIPTION
#   Here I will have generic methods to create synaptic weights between connected neurons.
#

import matplotlib.pyplot as plt
import numpy
import numpy.random

# I had a difficult time getting things to import, so now I have however, go about fixing this by going up a directoy
# which I think I can do in the following way.
import os
import sys
import pickle
sys.path.append(os.path.abspath("../"))

def uniformSampling(weightLow, weightHigh, connectionArray):
    """
    DESCRIPTION
        Here we will make synaptic weights by sampling over some uniform distribution.
    :param weightLow: The minimum value over the domian of the uniform distribution
    :param weightHigh: The maximum value over the domain of the unifrom distribution
    :param connectionArray: numpy.ndarray, nX4, in each row [[presynaptic Neuron ID, postsynaptic neuron ID, presynatpic neuron group,
        postsynaptic neuron group] ...]. This numpy.array is how we will store connections between neuron groups.
    :return: numpy.array(nX1) where n is equal to the number of rows in the given 'connectionArray'
    """
    weightArray = numpy.random.uniform(weightLow, weightHigh, (connectionArray.shape[0], 1))
    return weightArray