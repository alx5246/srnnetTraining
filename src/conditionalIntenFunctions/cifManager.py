# A. Lonsberry
# September 2016 - Taken from sandbox (condIntenFunc_0_3.py)
#
# DESCRIPTION
#   The class(es) below are designed to handle generation of a filter-manger class object, and a number of cifNeuron
#   objects.


import numpy
import matplotlib.pyplot as plt
import numpy.matlib
import time

import cifFilterManager
import cifNeuronCif

class condItenFuncManager:
    """
    DESCRIPTION
    The class is an overarching class that works a level above the classes condInenFuncFilterManager and conItenFunc.
    This is the class that sends data between the two of them and handles high level operations.

    GENERAL USE GUIDELINES:

    (1) Create new CIFs for neurons by calling self.addNewFilter(). This will internally handle creating the neuron's
        CIFs and handling all filter/function operations.

    (2) Initialize filter recorders by calling self.initFilterRecorders(), which will initialize all the necessary
        internal structures to save and store data.

    (3) Input new data from simulation using self.parseNewSpikeTrainData()

    (4) Prepare data for use in prediction of CIF updates using self.prepareParsedData()

    (5) ... call individual CIFs to update coefficients or make predictions .... this should be automated eventually

    """


    def __init__(self, dt):
        """
        self.filterManager: and instantiation of the condInenFuncFilterManger class
        self.condIntenFunctions: a list which will store all the conItenFunc class objects
        self.dt: double, the time-step we want to use for the CIFs
        """
        # An instantiattion of the filter-manager class
        self.filterManager = cifFilterManager.condInenFuncFilterManager()
        # Where we keep all our individual conditional intensity function objects for each individual neurons
        self.condIntenFunctions = []
        # The time-step we expect all the filters to be sampled at!!  I NEED TO REWRITE A BUNCH OF STUFF SO THAT WHEN
        # NEW FILTES ARE GIVEN, dt IS NO LOGER CONSIDERED OR NEEDED, RATHER WE SIMPLY REFERENCE THE GLOBAL VALUES HERE!
        self.dt = dt


    def getNeuronCifObj(self, neuronGroup, neuronId):
        """
        DESCRIPTION
        Determine if a particular neuron, denoted by its 'neuronGroup' and 'neuronId', has a CIF object within the list
        self.condIntenFunctions.

        :param neuronGroup: integer >= 0
        :param neuronId: interger >= 0

        :return: a list of objects, that may be empty if no CIF matches the given input
        """
        subList = [x for x in self.condIntenFunctions if x.neuronGroup==neuronGroup and x.neuronId==neuronId]
        return subList[0]


    def doesNeuronHaveCifObj(self, neuronGroup, neuronId):
        '''
        DESCRIPTION
        Determine if the given neuron (neuronId) from a particular neuron-group (neuronGroup) has a CIF class object
        instantiated

        :param neuronGroup: integer, the neuron-group the neuron is within
        :param neuronId: integer, the neuron ID

        :return: True -or- False
        '''
        subList = [x for x in self.condIntenFunctions if x.neuronGroup==neuronGroup and x.neuronId==neuronId]
        if len(subList)==0:
            toReturn = False
        else:
            toReturn = True
        return toReturn


    def createNeuronCifObj(self, neuronGroup, neuronId, dt):
        """
        DESCRIPTION
        Create a new CIF function object for some neuron.

        :param neuronGroup: integer, the neuron-group the neuron is within
        :param neuronId: integer, the neuron ID
        :param dt: double, time-step we want to use for the CIFs
        :return: the CIF object created, the last instantiation of conItenFunc in the self.condIntenFunctions list
        """
        self.condIntenFunctions.append(cifNeuronCif.condItenFunc(neuronGroup, neuronId, dt))
        return self.condIntenFunctions[-1] #Return the last created cifObj, which is the one we last created.


    def addNewFilter(self, newFilter):
        '''
        DESCRIPTION
        Adds a new function to neurons. More specifically, the input (see input description below)
        describes a type of filter/function for a CIF, for a neuron or neurons. The class here then alters or adusts
        its list of CIF objects (self.condIntenFunctions) and its filter-manager object (self.filterManager).

        We will outline the process,
        (1) Iterates over the given neurons in the 'newFilter' input
        (2) Checks to make sure the neuron has its own CIF object using the functions self.doesNeuronHaveCifObj,
            self.getNeuronCifObj, and self.createNeuronCifObj
        (3) Depending on the type of new filter/function, for example 'gausLinSyn', different things occur. In the most
            complex case where the CIF uses fitlered spike trains from other neuson, we iterate over all the synapses.
        (4) Check to make sure the new filter/function does not already
            eixist in the neurons CIF.
        (5) Call the filter manager objects add filter,
            self.filterManager.addFilter( ... ) which will (a) check to see if the filter exists, and (b) if it does not
            then it adds the new filter to the list of filter-types, and (c) add a new convlution filter if the filter
            is not simply a duplicate of a former filter with just a different time-delay.
        (6) Again call filter manager, and add the filter and neuron that needs to be filtered,
            self.filterManager.addNeuronAndFilter
        (7) Lastly update the CIF object to include the filter-type, the filter-postion that we need to reference in
            the self.filtereManger object, and the time-step delay.
        :param newFilter: is a list that can take the following forms,
               ['bias', neuronGroup, neuronIds] where neuronIds is a 1-D numpy.array object
               ['gausLinSelf', neuronGroup, neuronIds, a, b, c, timeLength, dt, timeDelays] where neuronIds and timeDelays are both 1-D numpy.array
               ['gausLinSyn', neuronGroup, neuronIds, synStruc, a, b, c, timeLength, dt, timeDelays]
               ['gausQuadSynSelf', neuronGroup, neuronIds, synStruc, a, b, c, timeLength, dt, timeDelays]
               NOTE, the acceptable form for CIF objects function-type is of the following forms,
               ['bias']
               ['gausLinSelf', a, b, c, timeLength, dt, timeDelay]
               ['gausLinSyn', otherNeuronGroup, otherNeuronId, a, b, c, timeLength, dt, timeDelay]
               ['gausQuadSynSelf', otherNeuronGroup, otherNeuronId, a, b, d, timeLength, dt, timeDelaySelf, timeDelayOther]

        TODO:
        (1) Add actual quadratic functions, they are currently not implemented

        :return: N/A
        '''
        # Iterate over the nerouns (neuronIds) in the given neuron-group (neuronGroup)
        for neuId in newFilter[2]:
            # Does a CIF object exist for the specifc neuronGroup + neuronId? If not, make a new CIF object for the neuron.
            if self.doesNeuronHaveCifObj(newFilter[1], neuId):
                cifObj = self.getNeuronCifObj(newFilter[1], neuId)    # Grab the object that already exists
            else:
                cifObj = self.createNeuronCifObj(newFilter[1], neuId, self.dt) # Create a new object if we need to
            # If the 'newFilter' input list includes the need to iterate over synaptic connections, we diverge here.
            if newFilter[0]=='gausLinSyn' or newFilter[0]=='gausQuadSynSelf':
                # Iterate over the 'synStruc'. The synapse strucuture is a numpy.array of size Nx4, where in each row there
                # exists information of the following format, [ [presynapticId, postsynapticId, preSynapticGroup,
                # postSynapticGroup], ... ].
                synStruc = newFilter[3]
                for i in range(synStruc.shape[0]):
                    # Grab the connecting neurons that connect the the particular neuronId + neuronGroup if possible
                    if synStruc[i,0]==neuId and synStruc[i,2]==newFilter[1]:
                        otherNeuronId    = synStruc[i,1]
                        otherNeuronGroup = synStruc[i,3]
                    elif synStruc[i,1]==neuId and synStruc[i,3]==newFilter[1]:
                        otherNeuronId    = synStruc[i,0]
                        otherNeuronGroup = synStruc[i,2]
                    else:
                        otherNeuronId    = -1
                        otherNeuronGroup = -1
                    # Iterate over the 'timeDelays' if we have a connection to our particular neuron from another
                    # synaptically connected neuron.
                    if otherNeuronId!=-1:
                        if newFilter[0]=='gausLinSyn':
                            for tDelay in newFilter[-1]:
                                # Check if this function/filter already exists in the cifObj, for each different type of
                                # filter type we have a different scheme. For the case of 'gausLinSyn' the input we are
                                # looking for follows the format, ['gausLinSyn', otherNeuronGroup, otherNeuronId, a, b, c,
                                # timeLength, dt, timeDelay]
                                dummyFunctionType = [newFilter[0], otherNeuronGroup, otherNeuronId, newFilter[4], newFilter[5],
                                                     newFilter[6], newFilter[7], newFilter[8], tDelay]
                                if not cifObj.doesFunctionExist(dummyFunctionType):
                                    # This particular filter for this neuorn does not exist. It does not mean the filter
                                    # does not already exist however. But we can simply call the same function to take care
                                    # of this problem, we simply input the following list ['gaussian', a, b, c, timeLength,
                                    # dt, timeDelay] into filterManager.addFilter() function, and it will return a list with
                                    # the values [truncation, filterIndex, timeStepsDelayed].
                                    outVals0 = self.filterManager.addFilter( ['gaussian', newFilter[4], newFilter[5],
                                                                             newFilter[6], newFilter[7], newFilter[8],
                                                                             tDelay] )
                                    # Now we have to indicate that this other neuorn specified by 'otherNeuronId' and
                                    # 'otherNeuronGroup' is going to be filtered, by only "adding" the combination, the
                                    # output below is the index we care about
                                    filterPosition = self.filterManager.addNeuronAndFilter(otherNeuronGroup, otherNeuronId, outVals0[1], outVals0[2], tDelay)
                                    # Ammend the 'cifObj' so we know which filter to point to and the number of time-steps
                                    # that are needed to be delayed
                                    cifObj.addFunction(dummyFunctionType, filterPosition, outVals0[2])
            # If the 'newFilter' is one of looking at self produced spikes
            elif newFilter[0]=='gausLinSelf':
                for tDelay in newFilter[-1]:
                    # Check if this function/filter alreayd exists in thecifObj... no duplicates allowed
                    #['gausLinSelf', neuronGroup, neuronIds, a, b, c, timeLength, dt, timeDelays]
                    #['gausLinSelf', a, b, c, timeLength, dt, timeDelay]
                    dummyFunctionType = [newFilter[0], newFilter[3], newFilter[4], newFilter[5], newFilter[6], newFilter[7], tDelay]
                    if not cifObj.doesFunctionExist(dummyFunctionType):
                        # This particular filter function for this neuron does not exist. It does not mean the filter
                        # does not alreayd exist however. But we can simply call the one function to take care of
                        # of this problem...
                        outVals0 = self.filterManager.addFilter( ['gaussian', newFilter[3], newFilter[4],
                                                                  newFilter[5], newFilter[6], newFilter[7], tDelay])
                        # Now we have to indicate that this neuron is going to be filtered, to do this we only have to
                        # call the "adding" function which will resolve any problems automatically (no duplicates)
                        filterPosition = self.filterManager.addNeuronAndFilter(newFilter[1], neuId, outVals0[1], outVals0[2], tDelay)
                        # Ammend teh cifObj
                        cifObj.addFunction(dummyFunctionType, filterPosition, outVals0[2])
            # Fi the 'newFilter' is a simply a 'bias' term.
            elif newFilter[0]=='bias':
                dummyFunctionType = [newFilter[0]] # ['bias']
                # Make sure we do not have any duplicates!
                if not cifObj.doesFunctionExist(dummyFunctionType):
                    # How do I handle this situation? The most simple is not to make any additions or alteratiosn to the
                    # filters (self.filterManger obj), but rather simply in the CIF object for the particular neuorn,
                    # make some simple addtions.
                    cifObj.addFunction(dummyFunctionType, -1, -1) # We give values of -1 for 'filterPostition' and 'stepsDelayed'


    def parseNewSpikeTrainData(self, listOfSpikeObjs, listOfNeuronGroups, startTime, endTime):
        """
        DESCRIPTION
        This method (1) pushes the input data the filter-manager class (self.filterManger) by calling
        self.filterManger.filterSpikeTrainData(), which will automatically filter the data, and store it in the correct
        form inside self.filterManager, and (2) indicate to the CIF objects that new data is available to to used to
        predict spikes and used to train and update CIR coefficients.

        :param listOfSpikeObjs: list of brian2.SpikeMonitor objects
        :param listOfNeuronGroups: a list corresponding to listOfSpikeObjs group-IDs, this tells us where each
            SpikeMonitor object is from
        :param startTime: the start-time from where we want to begin looking at or adding data
        :param endTime: the end-time from where the simulation ended or where we do not want to include data from. Note
            the end-time is not included as a time-step in any of the functions or methods.

        :return: N/A
        """
        # FIRST, give the data to the spike train-manager to let it parse through the data, and filter what it needs to.
        self.filterManager.filterSpikeTrainData(listOfSpikeObjs, listOfNeuronGroups, startTime, endTime)
        # SECONDLY, now the spikes of the neurons need to be reflected in each of the CIF objects. This is done by by
        # iterating over the conditional intensity functions objects stored in the the list, self.condIntenFunctions
        t0 = time.clock()
        for cifObj in self.condIntenFunctions:
            # Determine if and where our neurons are represented.
            if cifObj.neuronGroup in listOfNeuronGroups:
                spkObjInd = listOfNeuronGroups.index(cifObj.neuronGroup)
                spkObj = listOfSpikeObjs[spkObjInd]
                # 'spkObj' is a brian2.SpikeMonitor object. It may have spike times that precceedes 'startTime'. As such
                # only consider those times >= startTime and <= endTime, and lastly are the neuron we are looking for
                timeCorrectedInds = numpy.where((spkObj.t>=startTime) & (spkObj.t<=endTime) & (spkObj.i==cifObj.neuronId))
                if timeCorrectedInds[0].size>0:
                    spikeTimes = spkObj.t[timeCorrectedInds[0]]
                else:
                    spikeTimes = numpy.array([])
                # Now give the spike times to the cifObj, now both the self.filterManger object, and the CIF objects
                # all have the same updates having accessed the same data. Making sure to give the global time-step
                # self.dt.
                cifObj.appendSpikeArray(spikeTimes, startTime, endTime, self.dt)

        totalTime = time.clock() - t0
        outStr = "\nThe total time to store data in CIFs = " + str(totalTime)
        #print(outStr)


    def prepareParsedData(self):
        """
        DESCRIPTION
        This method individually iterates through the CIFs and prepares data to be used for prediction or coefficient
        updates.

        :return: N/A
        """
        for cifObj in self.condIntenFunctions:
            # Now that all the data is given to the cifObj, now we want to prepare the data for prediction
            cifObj.createFeasibleFilteredDataSet(self.filterManager.filteredArrays)


    def learnCifCoeffs(self, learningRate=.01, iterations=1, method='SGD'):
        """
        DESCRIPTION
        Method to update and train coefficients of the CIFs. This method operates under the following assumptions: (1)
        all the CIFs considered are stored here in this class within self.condIntenFunctions, (2) all the training data
        used is already given and stored in self.filterManager.filteredArrays and all feasible data there will be used.

        IMPORTANT, make sure to call self.prepareParsedData() first as that is needed to run this method.

        :param iterations: how many times do we want to iterate over the given data
        :param method: not used now, but in the future should be implemented if we should want to have different
            coefficient updating methods.
        :param learningRate: SGD learning rate

        :return: N?A
        """
        for cifObj in self.condIntenFunctions:
            cifObj.updateCifCoeffsSGD(learningRate)


    def garbageCollect(self, endTime):
        """
        DESCRIPTION
        This is a general garbage collection routine that is meant to be run after each time we (1)
        parseNewSpikeTrainData and (2) prepareParseData, and (3) learnCifCoeffs or predict spikes. Then this can be run
        which will dump all the data in the CIF objects and just hold onto the necessary data in the filter array's
        assuming the next data that will be given is from the next time steps in the simulation.
        :return:
        """
        self.filterManager.garbabgeCollectFilteredArrays()
        self.filterManager.garbageCollectSpikeTrainMemory(endTime)
        for cifObj in self.condIntenFunctions:
            cifObj.clearSpikeArray()


    def clearMemory(self):
        """
        DESCRIPTION
        Here ...
        :return:
        """
        self.filterManager.clearAllMemory()
        for cifObj in self.condIntenFunctions:
            cifObj.clearSpikeArray()


    def diagnosticCheckNeuronAndFilterCombs(self):
        """
        DESCRIPTION
        Here ...
        :return:
        """
        # (1) We want to look through all the conditional-instentiy objects stored
        print("\n")
        print("#######################################################################################################")
        print("#                                                                                                     #")
        print("#               Checking Conditional-Intensity Functions Stored in cif-manger obj                     #")
        print("#                                                                                                     #")
        print("#######################################################################################################")
        print(" ")
        if len(self.condIntenFunctions)<1:
            print("WARNING .... NO conditional-intensity function objects stored in manager obj")
        else:
            # How many CIF objects do we have?
            print("A total of",len(self.condIntenFunctions)," cif objects are stored in manager obj")
            print(" ")
            # Which neurons have CIF objects
            print("Displaying neuron and neuron-group combinations ...")
            for i in self.condIntenFunctions:
                print("... Neuron", i.neuronId, " of Neuron-Group", i.neuronGroup)
            # Look at what functions each CIF has
            print(" ")
            print("Looking at cifs for each neuron stored ... ")
            for cifObj in self.condIntenFunctions:
                print("... Neuron:", cifObj.neuronId, " of Neuron-Group:", cifObj.neuronGroup, " has the following cifs,")
                for fInd, funcType in enumerate(cifObj.functionType):
                    if len(funcType)>1: # Check if this a 'bias' type or not
                        outStr = "    ..."
                        for k in funcType:
                            outStr = outStr + " " + str(k) +","
                        print(outStr)
                        print("        ... which points to filtered spike train of neuron:",self.filterManager.neuronToFilter[cifObj.filterPositions[fInd],1],",of neuron group:",self.filterManager.neuronToFilter[cifObj.filterPositions[fInd],0],",in filter object")
                        print("        ... wherein the filter applied to this neuron has a time-step:",self.filterManager.filtersDt[self.filterManager.neuronToFilter[cifObj.filterPositions[fInd],2]])
                        print("        ... wherein the filter applied to this neuron has a backup-time:",self.filterManager.filtersBackUpTime[self.filterManager.neuronToFilter[cifObj.filterPositions[fInd], 2]])
                        outStr = "        ... whereing the filter applied to this has the function"
                        retFuncType = self.filterManager.getFilterType(self.filterManager.neuronToFilter[cifObj.filterPositions[fInd],2])
                        for k in retFuncType:
                            outStr = outStr + " " + str(k) + ","
                        print(outStr)
                        print("        ... where the number of time-stesp delayed is ",cifObj.filterStepsDelayed[fInd])
                    else:
                        outStr = "    ... " + funcType[0]
                        print(outStr)


    def diagnosticCheckExamineStoredFilters(self):
        print("\n")
        print("#######################################################################################################")
        print("#                                                                                                     #")
        print("#                                    Examining Store Filters                                          #")
        print("#                                                                                                     #")
        print("#######################################################################################################")
        print(" ")
        print("Looking at what filters are stored, and their associated filter vectors. The important point here is to\nmake sure that all the filters look correct")
        print(" ")
        for i, filtIndens in enumerate(self.filterManager.filterIdentities):
            plt.figure(i)
            if filtIndens[0]=='gaussian':
                outStr = "Convolution Filter " +  str(i) + "\na=" + str(filtIndens[1]) + ", b=" + str(filtIndens[2]) + ", c=" + str(filtIndens[3]) + ", t-len=" + str(filtIndens[4]) + ", trunc=" + str(filtIndens[7]) +  ", t-delay=" + str(filtIndens[6])
                plt.title(outStr)
            plt.plot(self.filterManager.filters[filtIndens[-1]])
        plt.show()


    def diagnosticCheckExamineFilteredSpikeData(self):
        print("\n")
        print("#######################################################################################################")
        print("#                                                                                                     #")
        print("#                                    Examining Filtered Data                                          #")
        print("#                                                                                                     #")
        print("#######################################################################################################")
        print(" ")
        print("Looking at what filters are stored, and their associated filter vectors. The important point here is to\nmake sure that all the filters look correct")
        print(" ")
        for i, filtArray in enumerate(self.filterManager.filteredArrays):
            plt.figure(i)
            outStr = "Filtered Spike Trian Data\n" + "Neuron ID=" + str(self.filterManager.neuronToFilter[i][1]) \
                     + ", Neuron Grp ID=" + str(self.filterManager.neuronToFilter[i][0]) + "\n Delay-Steps=" \
                     + str(self.filterManager.numbStepsToStore[i,0]) + ", Actual-Time-Delay=" + str(self.filterManager.actualDelay[i,0])
            plt.title(outStr)
            # Plot the filtered data, where data exists for each time-step!
            plt.plot(self.filterManager.filteredArraysTimeValues[i],filtArray)
            # Plot the actual spikes for the given neuron!
            times = numpy.matlib.repmat(self.filterManager.spikeTrainMemory[i],3,1)
            times = numpy.reshape(numpy.transpose(times),3*self.filterManager.spikeTrainMemory[i].size)
            vals = numpy.vstack((numpy.zeros(self.filterManager.spikeTrainMemory[i].size),numpy.ones(self.filterManager.spikeTrainMemory[i].size),numpy.zeros(self.filterManager.spikeTrainMemory[i].size)))
            vals = numpy.reshape(numpy.transpose(vals),3*self.filterManager.spikeTrainMemory[i].size)
            plt.plot(times, vals)
        plt.show()


    def diagnosticCheckExamineConstrutedDataSet(self):
        print("\n")
        print("#######################################################################################################")
        print("#                                                                                                     #")
        print("#                                    Examining Constructed Data                                       #")
        print("#                                                                                                     #")
        print("#######################################################################################################")
        print(" ")
        print("Looking at data compiled in the constructred data feild. As the data gets printed out to plots, we will\n"
              "also be plotting data out about the s")
        print(" ")
        for i, cifObj in enumerate(self.condIntenFunctions):
            plt.figure(i)
            outStr = "CIF Input Data\nNeuron=" + str(cifObj.neuronId) + ", Neuron-Grp=" + str(cifObj.neuronGroup)
            plt.title(outStr)
            for columnInd in range(cifObj.dataMatrix.shape[1]):
                # Note the data stored in the 'dataMatrix' numpy.array is such that the earliest time is given first, and
                # the lastest time is last, we need to flip the data so it shows up correctly
                plt.plot(cifObj.spikeTimeArray[-cifObj.dataMatrix.shape[0]:], cifObj.dataMatrix[:, columnInd], label=str(columnInd))
            plt.legend()
            plt.plot(cifObj.spikeTimeArray[-cifObj.dataMatrix.shape[0]:], cifObj.spikeArray, '--', linewidth=3.0)
        plt.show()

########################################################################################################################
#                                              Internal Unit Testing
########################################################################################################################


if __name__ == "__main__":


    ################# Testing the filter manager bits ##################################################################
    '''
    # Create a new conditional-intentisty-function manager object
    myFilterMan = condItenFuncManager(.001)
    # Create neurons and neuron groups
    neuronGroup0Ids = numpy.array([0,1,2,3])
    neuronGroup1Ids = numpy.array([0,1,2])
    synapses0       = numpy.array([[0,0,0,1],[0,1,0,1],[0,2,0,1],[1,2,0,1]]) # [ [nId1, nId2, nGrp1, nGrp2], .... ]
    gausFiltDelays  = numpy.array([0., .1, 1.])
    # Create new conditional intensity functiosn for neurons
    inputFunctionType = ['gausLinSyn', 1, neuronGroup1Ids, synapses0, 0., .5, .3, .1, .001, gausFiltDelays] #['gausLinSyn', neuronGroup, neuronIds, synStruc, a, b, c, timeLength, dt, timeDelays]
    myFilterMan.addNewFilter(inputFunctionType)
    #['gausLinSelf', neuronGroup, neuronIds, a, b, c, timeLength, dt, timeDelays] where neuronIds and timeDelays are both 1-D numpy.array
    inputFunctionType = ['gausLinSelf', 1, neuronGroup1Ids, 0., .5, .3, .1, .001, gausFiltDelays]
    myFilterMan.addNewFilter(inputFunctionType)
    # ['bias', neuronGroup, neuronIds] where neuronIds is a 1-D numpy.array object
    inputFunctionType = ['bias', 1, neuronGroup1Ids]
    myFilterMan.addNewFilter(inputFunctionType)
    # ['bias', neuronGroup, neuronIds] where neuronIds is a 1-D numpy.array object
    inputFunctionType = ['bias', 1, neuronGroup1Ids]
    myFilterMan.addNewFilter(inputFunctionType)
    myFilterMan.diagnosticCheckNeuronAndFilterCombs()
    '''


    ################# Testing the addition, and filtering of data ######################################################
    # Create a new conditional-intensity-function manager object
    myFilterMan = condItenFuncManager(.001) # We have to give it the dt we are using
    # Create neurons and neurons groups
    neuronGroup0Ids = numpy.array([0])
    neuronGroup1Ids = numpy.array([0,1])
    synapses0       = numpy.array([[0,1,0,1], [0,1,1,1]])
    #synapses0       = numpy.array([[0,1,1,1]])
    gausFiltDelays  = numpy.array([0., .005, .01, .015, .02])
    # Create new conditional intensity functison for neurons
    inputFunctionType = ['gausLinSyn', 1, neuronGroup1Ids, synapses0, 1., .0, .002, .1, .001, gausFiltDelays] #['gausLinSyn', neuronGroup, neuronIds, synStruc, a, b, c, timeLength, dt, timeDelays]
    myFilterMan.addNewFilter(inputFunctionType)
    # Create new conditional intensity functison for neurons
    inputFunctionType = ['bias', 1, neuronGroup1Ids]
    myFilterMan.addNewFilter(inputFunctionType)
    # Check the filters
    myFilterMan.diagnosticCheckNeuronAndFilterCombs()
    myFilterMan.diagnosticCheckExamineStoredFilters()
    #Now I create and add some data
    class someSpikeObj:
        def __init__(self, inds, times):
            self.t = times
            self.i = inds
    neuSpks = numpy.array([0, 0, 0, 0, 0, 0, 0])
    neuTmes = numpy.array([.001, .0081, .011, .025, .050, .0493, .08])
    neuTmes = neuTmes + .0001
    grp1SpikeMon = someSpikeObj(neuSpks, neuTmes)
    neuSpks = numpy.array([0, 1, 0, 1])
    neuTmes = numpy.array([.011, .025, .06, .075])
    neuTmes = neuTmes + .0001
    grp2SpikeMon = someSpikeObj(neuSpks, neuTmes)
    # Initialize the filter class
    myFilterMan.filterManager.initFilterRecorders()
    # Now actually give the data to be presented
    myFilterMan.parseNewSpikeTrainData([grp1SpikeMon, grp2SpikeMon], [0, 1], startTime=0., endTime=.125)
    # Now some debugging
    myFilterMan.diagnosticCheckExamineFilteredSpikeData()
    # Now add some new data and see what happends
    neuSpks = numpy.array([0, 0, 0, 0, 0, 0, 0])
    neuTmes = numpy.array([.001, .008, .011, .025, .050, .053, .08])
    neuTmes = neuTmes + .1251
    grp1SpikeMon = someSpikeObj(neuSpks, neuTmes)
    neuSpks = numpy.array([0, 1, 0, 1])
    neuTmes = numpy.array([.011, .025, .06, .075])
    neuTmes = neuTmes + .1251
    grp2SpikeMon = someSpikeObj(neuSpks, neuTmes)
    # Now actually give the data to be presented
    myFilterMan.parseNewSpikeTrainData([grp1SpikeMon, grp2SpikeMon], [0, 1], startTime=.125, endTime=.250)
    # Now some debugging
    #myFilterMan.diagnosticCheckExamineFilteredSpikeData()
    # Now we want to create the training data!
    myFilterMan.prepareParsedData()
    # Now look at all the training data in one plot
    myFilterMan.diagnosticCheckExamineConstrutedDataSet()
    # Now we want to try and train using gradient descent!
    def testSgdAndPrections(myFilterMan, trainIterations, cifObjNumb):
        # Grab the cif-object
        cifObj = myFilterMan.condIntenFunctions[cifObjNumb]
        # Begin to plot things by plotting out the original data!
        plt.figure(1)
        plt.plot(cifObj.spikeTimeArray[-cifObj.dataMatrix.shape[0]:], cifObj.spikeArray[-cifObj.dataMatrix.shape[0]:], '--', linewidth=3.0)
        # Generate initial predictions, and plot them
        probs = cifObj.generateSpikeProbabilities()
        plt.plot(cifObj.spikeTimeArray[-probs.size:], probs)
        #Iteritively train using SGD method
        for i in range(trainIterations):
            cifObj.updateCifCoeffsSGD(.01)
            probs = cifObj.generateSpikeProbabilities()
            plt.plot(cifObj.spikeTimeArray[-probs.size:], probs)
        plt.show()
    # Run the test method we just created
    testSgdAndPrections(myFilterMan, 500, 1)
