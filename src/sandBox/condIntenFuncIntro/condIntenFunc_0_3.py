# A. Lonsberry
# June 2017
#
# DESCRIPTION
#   This python file will have methods to create and train the conditional intensity functions that we use to predict
#   the spiking of the neurons.
#
#   In the simplest version we will have a single filter, we look at its value for different delays, this means we only
#   have to run a convolution operation for each neuron once!!! This should really speed up things.

import numpy
import matplotlib.pyplot as plt

########################################################################################################################
# Create Individual Elements of CIF
########################################################################################################################


def gaussianFilterGenerator(a, b, c, dt, startTime, endTime):
    '''
    DESCRIPTION
    Here we generate some Gaussian-like function discritized by the time-step 'dt'. More precisely, we create a gaussian
    like filter and then evaluate it at evenly spaced points and store this as a vector. This can then be used as one
    of the filters that we use to in some CIF.

    FUNCTIONAL FORM
    f(x) = a*exp(-(x-b)**2/(2*c)**2)

    :param a:
    :param b:
    :param c:
    :param startTime:
    :param stopTime:
    :param dt:
    :return:
    '''

    xvals = numpy.arange(start=startTime, stop=endTime, step=dt)
    convVector = a * numpy.exp(-1.*(((xvals-b)**2.)/(2.*c**2.)))
    return convVector


def halfGaussianFilterGenerator(a, b, c, dt, startTime, endTime):
    '''
    DESCRIPTION
    Here we generate some Gaussian-like function discritized by the time-step 'dt'. More precisely, we create a gaussian
    like filter and then evaluate it at evenly spaced points and store this as a vector. This can then be used as one
    of the filters that we use to in some CIF.

    The DIFFERNCE between this an your standard Gaussian, is that at input x>=0 we set the output to be = 0.

    '''

    xvals = numpy.arange(start=startTime, stop=endTime, step=dt)
    convVector = a * numpy.exp(-1. * (((xvals - b) ** 2.) / (2. * c ** 2.)))
    convVector[numpy.argwhere(xvals>=0)] = 0.
    return convVector

def truncatedGaussianFilterGenerator(a, b, c, dt, timeLength, truncation):

    xvals = numpy.arange(start=-1.*timeLength, stop=timeLength, step=dt)
    xvals = numpy.hstack((xvals, timeLength)) #Make sure we add the last time-step
    convVector = a * numpy.exp(-1. * (((xvals - b) ** 2.) / (2. * c ** 2.)))
    convVector[numpy.argwhere(xvals>=truncation)] = 0.
    return convVector


########################################################################################################################
# Convolutionl Functions
########################################################################################################################


def convolveSpikeTimes(convVector,spikeTimes,startTime,stopTime,dt):
    '''
    DESCRIPTION
    Given our spikeTime inputs, along with other bits, we calculate the values of the filter with the spike times.

    :param convVector:
    :param spikeTimes:
    :param startTime:
    :param stopTime:
    :param dt:
    :return:
    '''

    #timeArray = numpy.arange(start=startTime, stop=stopTime+dt, step=dt)
    boolArray = numpy.zeros(int((stopTime-startTime)/dt)+1) #We add one for the endpoint

    # Find indicies of the zeros array that we will want to set to a spike
    indsArray = ((spikeTimes - startTime)/dt)

    # Set the boolArray to have spikes where they should
    boolArray[indsArray.astype(int)] = 1.

    #Now convolve boolsArray with the convolution vector
    outputArray = numpy.convolve(boolArray, convVector, mode='full')

    #Output the entire array excpet the last values because they are not causel
    return outputArray[:-(convVector.size-1)]


########################################################################################################################
# Generation of Coefficients and Tie-Togethers
########################################################################################################################

# We need a set of algorithms to create and track what the coefficients we have are correlated too. MUCH of what will be
# here is perhaps overly hardcoded and needs to be broken up in future efforst.

class conItenFunc:
    """
    DESCRIPTION
    A class that contains CIF information for designated neurons that have such a CIF.

    self.neuronGroup:   integer, the group in which the neuron is from
    self.neuronId:      integer, the neuron index or ID associated to the neuron within some group
    self.functionType:  list of lists for the following forms,
                        ['bias']
                        ['gausLinSelf', a, b, c, timeLength, dt, timeDelay]
                        ['gausLinSyn', otherNeuronGroup, otherNeuronId, a, b, c, timeLength, dt, timeDelay]
                        ['gausQuadSynSelf', otherNeuronGroup, otherNeuronId, a, b, d, timeLength, dt, timeDelaySelf, timeDelayOther]
    self.fitlerPostions: list of indices, where each index points to a row of filtered spike data, wherein that data is
                         that which is needed to create the CIF
    self.filterStepsDelayed: the numeber of time-steps needed to look backwards in the filtered data vector
    """


    def __init__(self, neuronGroup, neuronId):
        self.neuronGroup        = neuronGroup
        self.neuronId           = neuronId
        self.functionType       = [] #List of function types, for example [ ["gausLinSelf', a, b, c, timeLength, dt, timeDelay], ...]
        self.filterPositions    = [] #The list of indicies that we grab the filtered data from, in order as the filters are created.
        self.filterStepsDelayed = [] #The list of delays used for each of the filtered data.


    def doesFunctionExist(self, newFunctionType):
        """
        DESCRIPTION
        Indicates if the given "newFunctionType" which can look like any of the lists acceptable for self.functionType,
        alreayd exists or does not exist
        :param newFunctionType: any of the acceptable lists for the self.functionType list.
        :return: True, False
        """
        dummyVal = False
        if len(self.functionType) > 0:
            for funcType in self.functionType:
                if funcType == newFunctionType:
                    dummyVal= True
                    break
        return dummyVal


    def getFunctionIndex(self, newFunctionType):
        """
        DESCRIPTION
        Find the index of the 'newFunctionType' in the self.functiionType variable
        :param newFunctionType: a list with of the same allowable format for self.functionType
        :return:
        """
        dummyVal = 0
        for i in range(len(self.functionType)):
            if self.functionType[i] == newFunctionType:
                dummyVal = i
                break
        return dummyVal


    def addFunction(self, functionType, filterPosition, stepsDelayed):
        """
        DESCRIPTION
        Add a new filter or functionType to the CIF object
        :param functionType:
        :param filterPosition:
        :param stepsDelayed:
        :return:
        """
        self.functionType.append(functionType)
        self.filterPositions.append(filterPosition)
        self.filterStepsDelayed.append(stepsDelayed)


class condInenFuncFilterManager:
    """
    DESCRIPTION
    Class: manages different filter-funtions and knows how to filter spike trains.
    """

    def __init__(self):
        """
        DESCRIPTION
        Conditional Intensity Filter Manager, this class manages all of the filters used to construct the CIFs for each
        given neuron. In this class there are 3 varaibles.

        self.filterIdentities: list of list, [ [filter 1 info], [filter 2 info], ... ],
                               where in each sublist the
                               following variants are acceptable,
                               ['gaussian', a, b, c, timeLength, dt, timeDelay, truncation, filterIndex]
        self.filters:          list of 1D, numpy.array, filters, where each filter is something we convolve the original
                               spike trians with.
        self.neuronToFilter    2D numpy.array, [ [neuron-group, neuron Id, filterIndex], ... ]
        """
        self.filterIdentities = []
        self.filters          = []
        self.neuronToFilter   = numpy.array([[]])


    def initFilterRecorders(self):
        """
        DESCRIPTION
        Once all the filters needed have been entered we need to initialize a list of numpy.arrays that will be the
        running filtered variables. It is only allowed to call this operator once.. for the time-being so we do not
        override earlier results...
        :return:
        """


    def getFilterIndex(self, filterInfo, timeDelay=None, truncation=None):
        """
        DESCRIPTION
        Given a filter input, it is determined if this ....
        :param filterInfo: a list of filter information, where the following varients are acceptable
                           ['gaussian', a, b, c, timeLength, dt]
        :return: filter index (integer) -or- None
        """
        toReturn = -1 #Initialize the output
        if filterInfo[0]=='gaussian':
            for aList in self.filterIdentities:
                if aList[0:6]==filterInfo:
                    if timeDelay!=None:
                        if aList[6]==timeDelay:
                            toReturn = aList[-1] #We always store the filterIndex last
                            break
                    else:
                        if aList[7]==truncation:
                            toReturn = aList[-1]  # We always store the filterIndex last
                            break
        return toReturn


    def getFilterTruncation(self, filterInfo, timeDelay=None, filterIndex=None):
        """
        DESCRIPTION
        Given a filter input, it is determined is this...
        :param filterInfo: a list of filter information, where the following varients are acceptable
                           ['gaussian', a, b, c, timeLength, dt]
        :return: filter truncation
        """
        toReturn = -1  # Initialize the output
        if filterInfo[0] == 'gaussian':
            for aList in self.filterIdentities:
                if aList[0:6] == filterInfo:
                    if timeDelay != None:
                        if aList[6] == timeDelay:
                            toReturn = aList[-2]  # We always store the filterIndex last
                            break
                    else:
                        if aList[8] == filterIndex:
                            toReturn = aList[-2]  # We always store the filterIndex last
                            break
        return toReturn


    def addFilter(self, filterInfo):
        """
        DESCRIPTION
        Add filters, an example being ['gaussian', a, b, c, timeLength, dt, timeDelay]. Now remember, there maybe more
        self.filterIdentities than actual self.fitlers. This is because many of those filters described in the
        self.filterIdentities are actually produce the same output just delayed by different time-steps!
        :param filterInfo: a list of filter information, where the following varients are acceptable,
                           ['gaussian', a, b, c, timeLength, dt, timeDelay]
        :return: [truncation, filterIndex]
        """
        if filterInfo[0]=='gaussian':
            #First determine the number of time-steps we need to look backwards.
            timeStepsDelayed = int(filterInfo[6]/filterInfo[5])
            if self.getFilterIndex(filterInfo[0:6], timeDelay=filterInfo[6]) != -1:
                # This filter already exists, which means we only have to grab the relevant output details which in this
                # case is the size of the filter's truncation, and the filter-index.
                dummyList = [self.getFilterTruncation(filterInfo[0:6],timeDelay=filterInfo[6]), self.getFilterIndex(filterInfo[0:6],timeDelay=filterInfo[6]), timeStepsDelayed]
            else:
                # Filter of that type does not exist, find truncation value and determine if we need to make an new
                # filter or if this is essentially another filter copy. Remember
                timeDelay  = filterInfo[6]
                timeLength = filterInfo[4]
                truncation = min(timeLength, timeDelay)
                if self.getFilterIndex(filterInfo[0:6], truncation=truncation) != -1:
                    #Filter already exists with different time-delay, (but really same filter)
                    newFilterIndex = self.getFilterIndex(filterInfo[0:6], truncation=truncation)
                    filterInfo.append(truncation)
                    filterInfo.append(newFilterIndex)
                    self.filterIdentities.append(filterInfo)
                else:
                    # In this case we have to create a new filter
                    self.filters.append(truncatedGaussianFilterGenerator(filterInfo[1], filterInfo[2], filterInfo[3], filterInfo[4], filterInfo[5], truncation))
                    newFilterIndex = len(self.filters)-1
                    filterInfo.append(truncation)
                    filterInfo.append(newFilterIndex)
                    self.filterIdentities.append(filterInfo)
                dummyList = [truncation, newFilterIndex, timeStepsDelayed]
        return dummyList


    def addNeuronAndFilter(self, neuronGroupId, neuronId, filterIndex):
        """
        DESCRIPTION
        Here we add
        :param neuronGroupId:
        :param neuronId:
        :param filterIndex:
        :return:
        """
        if self.neuronToFilter.size == 0:
            # It the original array is empty, then we need to make the first row
            self.neuronToFilter = numpy.array([[neuronGroupId, neuronId, filterIndex]])
            dummyOut = 0
        else:
            inds = numpy.where((self.neuronToFilter==(neuronGroupId, neuronId, filterIndex)).all(axis=1))
            if inds[0].size > 0:
                dummyOut = inds[0][0] #This is the position of the filtered values we are going to want to sample
            else:
                # The filter and neuron combo does not exist lets add it
                self.neuronToFilter = numpy.vstack(
                    (self.neuronToFilter, numpy.array([neuronGroupId, neuronId, filterIndex])))
                dummyOut = self.neuronToFilter.shape[0]-1 #Must subtract the 1 because we do zero indexing
        return dummyOut


class condItenFuncManager:
    """
    DESCRIPTION
    Class: manages all the conditional intensity functions for all the given neurons.
    """


    def __init__(self):
        """
        DESCRIPTION
        This is a class that manages all of the conditional intensity functions for all the neurons.
        """
        self.filterManager = condInenFuncFilterManager() # An instantiattion of the filter manager
        self.condIntenFunctions = [] # Where we keep all our individual conditional intensity function objects for each individual neuron


    def getNeuronCifObj(self, neuronGroup, neuronId):
        """
        DESCRIPTION
        Determine if a particular neuron, denoted by its 'neuronGroup' and 'neuronId', is represented in
        self.condIntenFunctions
        :param neuronGroup: integer >= 0
        :param neuronId: interger >= 0
        :return: the conItenFunc class instantiatition if it exists
        """
        subList = [x for x in self.condIntenFunctions if x.neuronGroup==neuronGroup and x.neuronId==neuronId]
        #Print returns a list... so I need to grab the first object of the list
        return subList[0]


    def doesNeuronHaveCifObj(self, neuronGroup, neuronId):
        '''
        DESCRIPTION
        Determine if the given neuron (neuronId) from a particular neuron-group (neuronGroup) has a CIF class object
        instantiated
        :param neuronGroup:
        :param neuronId:
        :return:
        '''
        subList = [x for x in self.condIntenFunctions if x.neuronGroup==neuronGroup and x.neuronId==neuronId]
        if len(subList)==0:
            toReturn = False
        else:
            toReturn = True
        return toReturn


    def createNeuronCifObj(self, neuronGroup, neuronId):
        self.condIntenFunctions.append(conItenFunc(neuronGroup, neuronId))
        return self.condIntenFunctions[-1] #Return the last created cifObj, which is the one we last created.


    def addNewFilter(self, newFilter):
        '''

        :param newFilter: is a list that can take the following forms,
               ['bias', neuronGroup, neuronIds] where neuronIds is a 1-D numpy.array object
               ['gausLinSelf', neuronGroup, neuronIds, a, b, c, timeLength, dt, timeDelays] where neuronIds and timeDelays are both 1-D numpy.array
               ['gausLinSyn', neuronGroup, neuronIds, synStruc, a, b, c, timeLength, dt, timeDelays]
               ['gausQuadSynSelf', neuronGroup, neuronIds, synStruc, a, b, c, timeLength, dt, timeDelays]


                ['bias']
                ['gausLinSelf', a, b, c, timeLength, dt, timeDelay]
                ['gausLinSyn', otherNeuronGroup, otherNeuronId, a, b, c, timeLength, dt, timeDelay]
                ['gausQuadSynSelf', otherNeuronGroup, otherNeuronId, a, b, d, timeLength, dt, timeDelaySelf, timeDelayOther]

        :return:
        '''
        # Iterate over the nerouns (neuronIds) in the given neuron-group (neuronGroup)
        for neuId in newFilter[2]:
            # Does a CIF object exist for the specifc neuronGroup + neuronId? If not, make one.
            if self.doesNeuronHaveCifObj(newFilter[1], neuId):
                cifObj = self.getNeuronCifObj(newFilter[1], neuId)    # Grab the object that already exists
            else:
                cifObj = self.createNeuronCifObj(newFilter[1], neuId) # Create a new object if we need to
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
                                filterPosition = self.filterManager.addNeuronAndFilter(otherNeuronGroup, otherNeuronId, outVals0[1])
                                # Ammend the 'cifObj' so we know which filter to point to.
                                cifObj.addFunction(dummyFunctionType, filterPosition, outVals0[2])


    def diagnosticCheckFilter(self, newFilter):
        """
        DESCRIPTION
        This is diagnostic, wherein here we run through some given filter to make sure all the components actually exist
        as they shoud.
        :param someFilter: is a list that can take the following forms,
               ['bias', neuronGroup, neuronIds] where neuronIds is a 1-D numpy.array object
               ['gausLinSelf', neuronGroup, neuronIds, a, b, c, timeLength, dt, timeDelays] where neuronIds and timeDelays are both 1-D numpy.array
               ['gausLinSyn', neuronGroup, neuronIds, synStruc, a, b, c, timeLength, dt, timeDelays]
               ['gausQuadSynSelf', neuronGroup, neuronIds, synStruc, a, b, c, timeLength, dt, timeDelays]
        :return:
        """
        # Iterate over the nerouns (neuronIds) in the given neuron-group (neuronGroup)
        for neuId in newFilter[2]:
            print("\n")
            print("Looking at CIF for neuron ",neuId," from neuron-group ",newFilter[1]," ... ")
            # Does a CIF object exist for the specifc neuronGroup + neuronId? If not, make one.
            if not self.doesNeuronHaveCifObj(newFilter[1], neuId):
                print("     WARNING : the neuron DOES NOT have a cifObj object stored")
            else:
                print("     the neuron Has a CIF object stored")
                cifObj = self.getNeuronCifObj(newFilter[1], neuId)
                # Iterate over the 'synStruc'. The synapse strucuture is a numpy.array of size Nx4, where in each row there
                # exists information of the following format, [ [presynapticId, postsynapticId, preSynapticGroup,
                # postSynapticGroup], ... ].
                synStruc = newFilter[3]
                for i in range(synStruc.shape[0]):
                    # Grab the connecting neurons that connect the the particular neuronId + neuronGroup if possible
                    if synStruc[i, 0] == neuId and synStruc[i, 2] == newFilter[1]:
                        otherNeuronId = synStruc[i, 1]
                        otherNeuronGroup = synStruc[i, 3]
                    elif synStruc[i, 1] == neuId and synStruc[i, 3] == newFilter[1]:
                        otherNeuronId = synStruc[i, 0]
                        otherNeuronGroup = synStruc[i, 2]
                    else:
                        otherNeuronId = -1
                        otherNeuronGroup = -1
                    # Iterate over the 'timeDelays' if we have a connection to our particular neuron from another
                    # synaptically connected neuron.
                    if otherNeuronId!=-1:
                        print("     looking at connection nId",synStruc[i,0],"nGrp",synStruc[i,2],",Into nId",synStruc[i,1],"nGrp",synStruc[i,3])
                        if newFilter[0] == 'gausLinSyn':
                            for tDelay in newFilter[-1]:
                                # Check if this function/filter already exists in the cifObj, for each different type of
                                # filter type we have a different scheme. For the case of 'gausLinSyn' the input we are
                                # looking for follows the format, ['gausLinSyn', otherNeuronGroup, otherNeuronId, a, b, c,
                                # timeLength, dt, timeDelay]
                                dummyFunctionType = [newFilter[0], otherNeuronGroup, otherNeuronId, newFilter[4],
                                                     newFilter[5],
                                                     newFilter[6], newFilter[7], newFilter[8], tDelay]
                                if not cifObj.doesFunctionExist(dummyFunctionType):
                                    print("          WARNING, the function filter does not exist!")
                                else:
                                    cifObjInd = cifObj.getFunctionIndex(dummyFunctionType)
                                    print("          funtion-type ",cifObj.functionType[cifObjInd][0],", filter position ",cifObj.filterPositions[cifObjInd],", steps delayed ",cifObj.filterStepsDelayed[cifObjInd])
        print("\n")

    def diagnosticAllFilters(self):
        if len(self.filterManager.filterIdentities)>0:
            for aFilter in self.filterManager.filterIdentities:
                print(aFilter)
        if len(self.filterManager.neuronToFilter)>0:
            for aFilter in self.filterManager.neuronToFilter:
                print(aFilter)


########################################################################################################################
# Internal Unit Testing
########################################################################################################################

if __name__ == "__main__":


    ##################Testing convolution filter creations

    #Test function generation
    #vectorVals = gaussianFilterGenerator(a=1.0, b=0., c=.01, dt=.001, startTime=-.1, endTime=.1)
    #plt.figure(1)
    #plt.plot(vectorVals)
    #vectorVals = halfGaussianFilterGenerator(a=1.0, b=0., c=.01, dt=.001, startTime=-.1, endTime=.1)
    #plt.plot(vectorVals)

    #Test the convolution to spike times thing
    #spikeTimes = numpy.array([1.001, 1.012, 1.045, 1.09, 1.231, 1.298, 1.521, 1.529, 1.545])
    #evaled = convolveSpikeTimes(vectorVals, spikeTimes, startTime=1.000, stopTime=1.6, dt=.001)
    #startTime = 1.6 - (evaled.size-1)*.001
    #print(startTime)
    #timesArray = numpy.arange(start=startTime, stop=1.6+.001, step=.001)
    #print(timesArray.shape)
    #plt.figure(2)
    #plt.plot(timesArray, evaled)
    #plt.show()


    #################Testing the filter manager bits
    myFilterMan = condItenFuncManager()
    #Create some neurons and neuron groups
    neuronGroup0Ids = numpy.array([0,1,2,3])
    neuronGroup1Ids = numpy.array([0,1,2])
    synapses0 = numpy.array([[0,0,0,1],[0,1,0,1],[0,2,0,1],[1,2,0,1]])
    timeDelays = numpy.array([0.,.01,.02,1.5,2.0])
    inputFunctionType = ['gausLinSyn', 1, neuronGroup1Ids, synapses0, 0, .5, .3, .1, .001, timeDelays]
    myFilterMan.addNewFilter(inputFunctionType)
    myFilterMan.diagnosticCheckFilter(inputFunctionType)
    myFilterMan.diagnosticAllFilters()

