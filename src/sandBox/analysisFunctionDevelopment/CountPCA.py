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

def count_PCA(countArray):
    """
    FUNCTION DESCRIPTION
        This function takes as input a n-D numpy.array (where n is the number of neurons in network) is the  of spike
        counts for each neuron, and outputs the principal component analysis (PCA) analyzed variance contribution


        This is an example function, and thus has some limitations. For example it can only handle 1D input. A more
        complex version of this function would likely be able to handle N-D input (spike trains from multiple neurons)
        simultaneously and would average the number of spikes over a user definesd interval and then average over the
        number of neurons so that each averaging interval, a scalar mean would be outputed. This function is not
        concerned with identifying how many spikes per interval are occurring for each neuron, but rather an averaged
        count across all neurons (i.e. the network) per interval

        Steps of PCA:


    :param countArray: n-D numpy.array, units are seconds, neuron spike times for each neuron stored in an numpy.array

    """

#Toy data for spike times (rows) for each neuron (columns)
# (i.e. row one has spike times for neuron one across the whole run)
#This toy day emulates spike times that might occur for a three neuron network
# countArray = ([.1, 1.3, 1.4, 2.2, 2.7, 3.5, 4.3, 4.4], [.1, 1.3, 1.4, 2.2, 2.7, 3.5, 4.3, 4.4], [.1, 1.3, 1.4, 2.2, 2.7, 3.5, 4.3, 4.4])
countArray = ([.1, 1.3, 1.4, 2.2, 2.4, 3.5, 4.4, 4.4], [0, 1.3, 1.3, 2.4, 2.9, 3.4, 4.3, 4.8], [.5, 1.8, 1.9, 2.3, 2.6, 3.5, 4.3, 4.6])


# Spike time turned into a numpy array (to ensure type)
countArray = np.array(countArray)
print(countArray)

# #One line of code?
results = PCA(countArray.T)
print ('Proportion of Variance: ', results.fracs, '\n')
print ('Eigenvalues: ', results.s, '\n')
print ('Weights: ', results.Wt, '\n')
Wt = np.array(results.Wt)
print ('Projected PCA: ', results.Y)
