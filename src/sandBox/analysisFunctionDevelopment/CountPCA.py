# S. Pickard
# June 2016
# (Using python 3.4.4, and brian2 2.0rc)
#
# DESCRIPTION
#   Function/method that produces an output
#   This analysis was inspired by Lazar et al.(2009) Figure 3.H
#   Goal is

import numpy as np
from matplotlib.mlab import PCA

def count_PCA(spikeTimeArray):
    """
    FUNCTION DESCRIPTION
        This function takes as input an N-D numpy.array (where n is the number of neurons in network) is the  of spike
        times for each neuron, and outputs the principal component analysis (PCA) analyzed variance contribution

    :param countArray: n-D numpy.array, units are seconds, neuron spike times for each neuron stored in an numpy.array

    """

#Toy data for spike times (rows) for each neuron (columns)
# (i.e. row one has spike times for neuron one across the whole run)
#This toy day emulates spike times that might occur for a three neuron network
# spikeTimeArray = ([.1, 1.3, 1.4, 2.2, 2.7, 3.5, 4.3, 4.4], [.1, 1.3, 1.4, 2.2, 2.7, 3.5, 4.3, 4.4], [.1, 1.3, 1.4, 2.2, 2.7, 3.5, 4.3, 4.4])
spikeTimeArray = ([.1, 1.3, 1.4, 2.2, 2.4, 3.5, 4.4, 4.4], [0, 1.3, 1.3, 2.4, 2.9, 3.4, 4.3, 4.8], [.5, 1.8, 1.9, 2.3, 2.6, 3.5, 4.3, 4.6])


# Spike time turned into a numpy array (to ensure type)
spikeTimeArray = np.array(spikeTimeArray)
print(spikeTimeArray)

# #One line of code?
results = PCA(spikeTimeArray.T)
print ('Proportion of Variance: ', results.fracs)


# print ('Eigenvalues: ', results.s, '\n')
# print ('Weights: ', results.Wt, '\n')
# Wt = np.array(results.Wt)
# print ('Projected PCA: ', results.Y)
