# S. Pickard
# June 2016
# (Using python 3.4.4, and brian2 2.0rc)
#
# DESCRIPTION


import numpy as np
from spikeMonToMatrix import spikeMon_To_Matrix
from SpikeCount import spike_count
import pickle
import matplotlib.pyplot as plt

########################################################################################################################
# TESTING PCA of spike counts in user defined sub time intervals

#Test has a series of functions that need to called prior to using the countPCA function. The following list is in
#sequential order of operations with a short description of what each function is doing:

#   1) spikeMonToMatrix: takes spikemon.i and spikmon.t arrays from simulation and generates an output array where the
#      rows correspond to neuron index and the columns are spike times
#   2) SpikeCount: counts the number of spikes in a given interval. Output is an array sized N x Int (number of neurons
#      by the number of sub time intervals).
#   3) countPCA: take in array and performs PCA on array of spike counts for each neuron
########################################################################################################################


# Get spikemon data for spike times from simulation and save to variable
inputFile = open("savedData_0/netOutput0_PoiNeu_SpikesTimes.pkl","rb")
spikeTimes = pickle.load(inputFile)
spikeTimesUnits = pickle.load(inputFile)
inputFile.close()
# print (spikeTimes)
type (spikeTimes)

# Get spikemon data for neuron indices from simulation and save to variable
inputFile = open("savedData_0/netOutput0_PoiNeu_SpikesInds.pkl","rb")
spikeTimeInds = pickle.load(inputFile)
inputFile.close()
# print (spikeTimeInds)


#With neuron indices and spike times array, generate matrix that holds the spike times of each neuron in each element
NeurFire = spikeMon_To_Matrix(spikeTimeArray = spikeTimes, NeurIndexArray = spikeTimeInds)
NeurFire = np.array(NeurFire)


# # Find the spike times for a particular neuron
# neuronSpikeTimes = spikeTimes[spikeTimeInds==99]
# # Now find the firing rate
# smoothedRate = instant_firing_rate(spikeTrain=neuronSpikeTimes, startTime=0.0, endTime=10.0, filterLength=5.0, var=1.)
# plt.figure(2)
# plt.plot(smoothedRate[1],smoothedRate[0])
# plt.show()

