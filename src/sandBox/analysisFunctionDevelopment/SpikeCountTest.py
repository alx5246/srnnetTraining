# S. Pickard
# June 2016
# (Using python 3.4.4, and brian2 2.0rc)
#
# DESCRIPTION
#   Here test a function/method that produces a smooth firing rate approximation of neuron spike output.
#
#   Here we will show how methods for looking at the firing rate of a single neuron over time. These methods are
#   enumerated in firingRateDistrib_0_0.py, but are called and tested in this file.
#
#   Using data that has been saved (perhaps saved in "savedData_0/" and generated by "createData_0_0.py") here we will
#   demonstrate some ways to look at firing rate statistics.


import numpy as np
from SpikeCount import spike_count
import pickle
import matplotlib.pyplot as plt


# import plotly.plotly as py
# import plotly.graph_objs as go

########################################################################################################################
# TESTING Spike Counting
########################################################################################################################
#This segment works
# In the first test we will simply create some fake neuron output a specific frequency.
# staticFiringPattern = np.linspace(start=0, stop=5, num=10)
# print (staticFiringPattern)
# spikeCount = spike_count(spikeTime=staticFiringPattern, start=0 , stop=5, dt=.3)
# print (spikeCount)

# In the second test we will load some data that we have created in simulation and see what it looks like
# Input spike times
inputFile = open("savedData_0/netOutput0_PoiNeu_SpikesTimes.pkl","rb")
spikeTimes = pickle.load(inputFile)
# spikeTimesUnits = pickle.load(inputFile)
inputFile.close()
# print (spikeTimes)

start = 0.0
stop = 10.0
dt = .2

# # Input spike times indices,
# inputFile = open("savedData_0/netOutput0_PoiNeu_SpikesInds.pkl","rb")
# # spikeTimeInds = pickle.load(inputFile)
# inputFile.close()
# # Find the spike times for a particular neuron
# neuronSpikeTimes = spikeTimes[spikeTimeInds==99]
#
# # Now find spike count
[spikeCount, time] = spike_count(spikeTime=spikeTimes, start=start, stop=stop, dt = dt)
spikeCount = np.array(spikeCount)
# print (spikeCount)
# print (time)

nSpike = len(spikeCount)
nTime = time_len = (len(time))
print (nTime)

bin_size = 0.2; min_edge = 0; max_edge = 10
N = (max_edge-min_edge)/bin_size; Nplus1 = N + 1
bin_list = np.linspace(min_edge, max_edge, Nplus1)
print (bin_list)

# plt.figure(1)
plt.hist(list(spikeCount), color="#3F5D7D", bins=bin_list)
# plt.bar(left, height = spikeCount, width = 1, facecolor = 'blue')
plt.title("Spike Count Across Time Intervals")
plt.xlabel("Time Intervals")
plt.ylabel("Spike Count")
# plt.axis([0,50])
plt.grid(True)
plt.show()
