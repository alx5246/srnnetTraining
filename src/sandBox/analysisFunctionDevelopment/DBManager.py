import pickle
import json
import os

class DBManager(object):

    #
    def ___init__(self, fileName):
        ##description
        self.fileName = fileName

        if os.path.exists(fileName) == True:
            with open(fileName, 'r') as file:
                self.fileData = json.load(file)
        else:
            self.fileData = dict()


    def record_Data(self, alpha, beta, tau, gamma, time, timestep, outputData):
        paramKey = "alpha=%s,beta=%s,tau=%s,gamma=%s" % (alpha, beta, tau, gamma)
        if paramKey not in self.fileData:
            self.fileData[paramKey] = {}
        # at this point, we have a dictionary reference to the parameters of the simulation

        timeKey = "time=[0, %s],timestep=%s" % (time, timestep)
        if timeKey not in self.fileData[paramKey]:
            self.fileData[paramKey][timeKey] = {}
        # at this point we have a dictionary reference to the time values of the simulation
        # FOR given parameters

        # outputData is a dictionary of all the pickled data run in the simulation keyed by the simulation
        # variables (K, R_hat, etc.)
        for pair in outputData.items(): # dict.items() returns a list of [key, value] pairs
            if pair[0] not in self.fileData[paramKey][timeKey]:
                self.fileData[paramKey][timeKey][pair[0]] = [pair[1]]
            else:
                self.fileData[paramKey][timeKey][pair[0]].append(pair[1])


    def save_Data(self):
        with open(self.fileName, 'w') as file:
            json.dump(obj = self.fileData, fp = file)


    def data_Dig(self, alpha, beta, tau, gamma, time, timestep):
        paramKey = "alpha=%s,beta=%s,tau=%s,gamma=%s" % (alpha, beta, tau, gamma)
        timeKey = "time=[0, %s],timestep=%s" % (time, timestep)

        outputDict = dict()
        for simVarPair in self.fileData[paramKey][timeKey].items(): # simVarPair = [k, R_hat, etc. : ["data", "units"]
            outputDict[simVarPair[0]] = dict()
            for dataPair in simVarPair[1].items(): # unpickle "data" entry and "units" entry
                outputDict[simVarPair[0]][dataPair[0]] = pickle.loads(dataPair[1])
        return outputDict