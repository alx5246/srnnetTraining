# S. Pickard
# June 2016
# (Using python 3.4.4, and brian2 2.0rc)
#
# DESCRIPTION


import numpy as np
from spikeMonToMatrix import spikeMon_To_Matrix
from SpikeCount2D import spike_count2D
from CountPCA import count_PCA
import pickle
import matplotlib.pyplot as plt

########################################################################################################################
# TESTING PCA of spike counts in user defined sub time intervals

#Test has a series of functions that need to called prior to using the countPCA function. The following list is in
#sequential order of operations with a short description of what each function is doing:

#   1) spikeMonToMatrix: takes spikemon.i and spikmon.t arrays from simulation and generates an output array where the
#      rows correspond to neuron index and the columns are spike times
#   2) SpikeCount2D: counts the number of spikes in a given interval. Output is an array sized N x Int (number of neurons
#      by the number of sub time intervals).
#   3) countPCA: take in array and performs PCA on array of spike counts for each neuron
########################################################################################################################




##################
#FAKE DATA
##################
# In the first test we will simply create some fake neuron output a specific frequency.

# staticFiringPattern = np.linspace(start=0, stop=5, num=10)
# staticFiringPattern = np.array(staticFiringPattern)
# staticFiringPattern = ([ 0., 0.55555556, 1.11111111, 1.66666667, 2.22222222, 2.77777778, 3.33333333, 3.88888889, 4.44444444, 5.], [ 0., 0.3, 1.5, 1.66666667, 2.8, 2.77777778, 3.5, 3.88888889, 4.6, 5.], [ 0., 0.3, 1.6, 1.66666667, 2.5, 2.77777778, 3.2, 3.9, 4.3, 5.])
# print (staticFiringPattern)


##################
#SIMULATED DATA
##################
# Get spikemon data for spike times from simulation and save to variable
inputFile = open("savedData_0/netOutput0_PoiNeu_SpikesTimes.pkl","rb")
spikeTimes = pickle.load(inputFile)
spikeTimesUnits = pickle.load(inputFile)
inputFile.close()
# print (spikeTimes[0,:])
# print (spikeTimes[1])
# print (spikeTimes[2])
# print (len(spikeTimes))

# Get spikemon data for neuron indices from simulation and save to variable
inputFile = open("savedData_0/netOutput0_PoiNeu_SpikesInds.pkl","rb")
spikeTimeInds = pickle.load(inputFile)
inputFile.close()
# print (len(spikeTimeInds))
# print (spikeTimeInds)
# print (max(spikeTimeInds))

#With neuron indices and spike times array, generate matrix that holds the spike times of each neuron in each element
NeurFire = spikeMon_To_Matrix(spikeTimeArray = spikeTimes, NeurIndexArray = spikeTimeInds)
NeurFire = np.array(NeurFire)
# print (NeurFire)
# print(len(NeurFire))


##################
#FIRE COUNT MATRIX
##################
# Now find spike count
start = 0.0
stop = 10.0
dt = 1

[spikeCount, timeInt] = spike_count2D(spikeTime=NeurFire, start=start , stop=stop, dt=dt)
spikeCount = np.array(spikeCount)
print (spikeCount)
# print (len(spikeCount))
# spikeCount.shape
print (timeInt)

##################
#PCA ON FIRE COUNT
##################

PCAcount = count_PCA(spikeCountArray = spikeCount)
print ('Proportion of Variance: ', PCAcount.fracs*100)
# print('Eigenvalues: ', PCAcount.s, '\n')
# print('Weights: ', PCAcount.Wt, '\n')

#Plot Principal components in histogram
a = PCAcount.fracs

data = range(10)
data = np.array(data)

d = np.diff(np.unique(data)).min()
left_of_first_bin = data.min() - float(d)/2
right_of_last_bin = data.max() + float(d)/2
interval = [1,2,3,4,5,6,7,8,9,10]
cum_var_exp = np.cumsum(a)
plt.hist(interval, weights = a)
plt.step(range(1,11), cum_var_exp, where='mid',
         label='cumulative explained variance')
axes = plt.gca()
axes.set_xlim(1, 10)
axes.set_ylim(0, 1)
axes.set_xticklabels([.1])
axes.set_yticklabels([])
plt.show()






major_ticks = np.arange(0, 1, .1)

plt.stem(interval,a)
plt.step(range(1,11), cum_var_exp, where='mid', label='cumulative explained variance')
axes = plt.gca()
axes.set_xlim(1, 10)
axes.set_ylim(0, 1)
axes.set_xticklabels([.1])
axes.set_yticklabels([1])
plt.show()