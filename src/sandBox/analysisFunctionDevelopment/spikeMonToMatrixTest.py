# S. Pickard
# June 2016
# (Using python 3.4.4, and brian2 2.0rc)
#
# DESCRIPTION



import numpy as np
from spikeMonToMatrix import spikeMon_To_Matrix
import pickle
import matplotlib.pyplot as plt

import os
import sys


currentDir = os.path.dirname(os.path.realpath(__file__)) # returns a string
layerStripping = [".." for x in range(0, 4)] # a list of ".." strings ["..", "..", "..", ".."]

# RELATIVE TO THE DIRECTORY controller.py is in...currentDir is the string form
# of the file path to the directory that controller.py is in...modify it
dbPath = os.path.join(currentDir, *layerStripping) # equivalent of os.path.join(currentDir, "..", "..", "..", "..")

# let python know where to search
sys.path.extend([currentDir, dbPath])

import db

manager = db.DBManager(databaseName="carlsonOneNeuron")
manager.openCollection("STDParameterRange")

########################################################################################################################
# TESTING spike monitoring data to matrix
########################################################################################################################

def getDataOut(recordList):
    recordDataList = []
    recordUnitList = []
    # get other information back like alpha, beta, gamma, etc.
    # remember you have access to ALL the data in the records from the database!
    for record in recordList:
        recordDataList.append(pickle.loads(record['data']['data']))
        recordUnitList.append(pickle.loads(record['data']['units']))
    return recordDataList, recordUnitList

# # In the second test we will load some data that we have created in simulation and see what it looks like
# # Input spike times
# inputFile = open("savedData_0/netOutput0_PoiNeu_SpikesTimes.pkl","rb")
spikeTimes, spikeTimesUnits = manager.query(
    {'measuredVariable': 'Neu_t'},
    hook=getDataOut
)
# spikeTimes = pickle.load(inputFile)
# spikeTimesUnits = pickle.load(inputFile)
# inputFile.close()
# print (spikeTimes)
# print (len(spikeTimes))

# Input spike times indices,
# inputFile = open("savedData_0/netOutput0_PoiNeu_SpikesInds.pkl","rb")
spikeTimeInds = manager.query(
    {'measuredVariable': 'Neu_i'},
    hook=getDataOut
)[0]
# print (spikeTimeInds)
# print (len(spikeTimeInds))


alpha = manager.query(
    {'paramters': 'alpha'},
    hook=getDataOut
)

# print (alpha)



assert(len(spikeTimeInds) == len(spikeTimes))
assert(len(spikeTimeInds) == len(spikeTimesUnits))


for simulationIndex in range(0, len(spikeTimeInds)):
    NeurFire = spikeMon_To_Matrix(spikeTimeArray = spikeTimes[simulationIndex],
                                  NeurIndexArray = spikeTimeInds[simulationIndex])
    NeurFire = np.array(NeurFire)
    # nSpike = len(spikeCount)
    # nTime = time_len = (len(time))
    print (NeurFire)








# nump_hist = plt.figure()
# plt.hist(spikeCount, bins = time)
# plt.figure(1)
# plt.hist(spikeCount, bins=50)
# # plt.bar(left, height = spikeCount, width = 1, facecolor = 'blue')
# plt.title("Spike Count Across Time Intervals")
# plt.xlabel("Time Intervals")
# plt.ylabel("Spike Count")
# # plt.axis([0,50])
# plt.grid(True)
# plt.show()

