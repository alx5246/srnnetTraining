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
import numpy.matlib

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

def truncatedGaussianFilterGenerator(a, b, c, timeLength, dt, truncation):
    """
    DESCRIPTION
    Create a truncated Gaussian-type filter. Things are truncated so when recorded spike trains are filtered, spikes
    forward in time do not get added to the filter as necessary, and thus keeping data points causel.

    Input 'truncation' is measured in seconds and indicates the first piont in time where the data should be recorded.
    Thus if 'timeLength' is 1.0 seconds which means the Gaussian filter spans from -1 to 1 seconds, and the truncation
    is 0.0, this means all the points in the fitler vector >0 are made to 0.

    IMPORATANT TO NOTE, that because we have a truncation, we are inplicitly storing a time-delay in the filter!

    :param a:
    :param b:
    :param c:
    :param timeLength:
    :param dt:
    :param truncation:
    :return:
    """

    xvals = numpy.arange(start=-1.*timeLength, stop=timeLength+dt, step=dt)  #Make sure we add the last time-step by adding dt to 'stop'
    convVector = a * numpy.exp(-1. * (((xvals - b) ** 2.) / (2. * c ** 2.)))
    # We want the value at the time = truncation to have a value, but everything greater has to be equal to zero. And
    # more specifically we just get rid of those extra values by deleting them from the array
    convVector = convVector[xvals<truncation]
    # When later calling numpy.convole, the vector get flipped along the axis, so we pre-flip here so when called it
    # gets flipped back as expected
    convVector = convVector[::-1]
    #print('hey')
    return convVector


########################################################################################################################
# Convolutionl Functions
########################################################################################################################


def convolveSpikeTimes(convVector, spikeTimes, startTime, stopTime, dt):
    '''
    DESCRIPTION
    Given our spikeTimes, we need to to convle the convolution filter with the the spikes. This is done by making a
    boolean array and then convolving with that array. IMPORTANT to note, we do NOT include the timing overlaps,
    that is we return a filtered vector from 'startTime' to the given 'stopTime'. IMPORTANT not note, we include the
    given stopTime as well.

    :param convVector:
    :param spikeTimes:
    :param startTime:
    :param stopTime:
    :param dt:
    :return:
    '''

    timeArray = numpy.arange(start=startTime, stop=stopTime+dt, step=dt) #This is going to be for the output
    boolArray = numpy.zeros(int((stopTime-startTime)/dt)+1) #We +1 in order to add the end time.
    # Find indicies of the zeros array that we will want to set to a spike
    indsArray = ((spikeTimes - startTime)/dt)
    # Set the boolArray to have spikes where they should
    boolArray[indsArray.astype(int)] = 1.
    #Now convolve boolsArray with the convolution vector
    outputArray = numpy.convolve(boolArray, convVector, mode='full')
    #Now trimp the array down both side of the data
    finalOutputArray = outputArray[(convVector.size-1):-(convVector.size-1)]
    #Output the entire array excpet the last values because they are not causel
    return [finalOutputArray, timeArray]

def convolveSpikeTimesWthBackUp(convVector, spikeTimes, startTime, stopTime, dt, backupTime=0):
    """
    DESCRIPTION
    Convolve a boolean spike array from either (1) startTime to stopTime or (2) startTime-backupTime to stopTime. In
    either case, the output vector is only going to be over a time interval from [startTime, stopTime]. The inclusion of
    backupTime, is to make the convolution more accurate at the edges.

    :param convVector:
    :param spikeTimes:
    :param startTime:
    :param stopTime:
    :param dt:
    :param backupTime:
    :return: [filteredArray, timeArray] where filteredArray is a vector of filtered results at each time-step that
             over the the doamian [startTime, stopTime], and timeArray is also a 1D numpy.array of equal length but
             each part of the vector includes the corresponding time-step value.
    """
    # Use numpy.arange and do not include the last time-step
    timeArray = numpy.arange(start=startTime-backupTime, stop=stopTime, step=dt)  # This is going to be for the output
    boolArray = numpy.zeros(int((stopTime - (startTime-backupTime)) / dt))
    indsArray = ((spikeTimes - (startTime-backupTime))/dt)
    # Make sure if there is a spike that has a spike-time of end-time, we make it slightly less so it rounds down
    indsArray[indsArray>=boolArray.size] -= 1
    boolArray[indsArray.astype(int)] = 1.
    outputArray = numpy.convolve(boolArray, convVector, mode='valid')
    numbStepsToGrab = int((stopTime-startTime)/dt)
    finalOutputArray = outputArray[-numbStepsToGrab:]
    timeArray = timeArray[-numbStepsToGrab:]
    return [finalOutputArray, timeArray]

########################################################################################################################
# Generation of Coefficients and Tie-Togethers
########################################################################################################################

# We need a set of algorithms to create and track what the coefficients we have are correlated too. MUCH of what will be
# here is perhaps overly hardcoded and needs to be broken up in future efforst.

class conItenFunc:
    """
    DESCRIPTION
    A class that contains CIF information for designated neurons that have such a CIF.

    General Guidelines on useage:
    (1) Add fucntions to the CIF object: typtically after an object of this class is created, you can add the functions
        that comprise the GLM suing the self.addFunction() method. A CIF as it is here described is a GLM that is passed
        through a logistic function. Everytime one uses self.addFunction() mehtod you are adding another linear component
        to the GLM
    (2) Add spikes and create data: once a simulation has been run and data is available there are two methods that need
        to be called. First is the self.appendSpikeArray() which takes in the new spikes and adds them to the memory
        here. Second is the self.createFeasibleFitleredDataSet() which takes in data from  a filter-manager object and
        creates a matrix of data that can be used to generate probabilities.
    (3) Generate ...
    """


    def __init__(self, neuronGroup, neuronId, dt):
        """
        self.neuronGroup: integer, the group in which the neuron is from
        self.neuronId: integer, the neuron index or ID associated to the neuron within some group
        self.functionType:  list of lists for the following forms,
                            ['bias']
                            ['gausLinSelf', a, b, c, timeLength, dt, timeDelay]
                            ['gausLinSyn', otherNeuronGroup, otherNeuronId, a, b, c, timeLength, dt, timeDelay]
                            ['gausQuadSynSelf', otherNeuronGroup, otherNeuronId, a, b, d, timeLength, dt, timeDelaySelf, timeDelayOther]
        self.filterPostions: list of indices, where each index points to a row of filtered spike data, wherein that data is
                             that which is needed to create the CIF. NOTE, if this value is -1, this means we do NOT search
        self.filterStepsDelayed: the numeber of time-steps needed to look backwards in the filtered data vector
        self.filterCoefficients: 2D numpy.array, Nx1, where N is the number of functions.
        self.dataMatrix: nXm 2D numpy.array, where the number of rows is equal to the number of data samples, and the
                         number of columns equal to the number of coefficients. When this is created, it is ASSUMED,
                         that the data samples come from consecutive time-steps, with the last data-point coming from
                         the point in time last recorded here in self.spikeBoolTimeArray.
        self.spikeBoolArray:

        :param neuronGroup:
        :param neuronId:
        """
        self.neuronGroup        = neuronGroup
        self.neuronId           = neuronId
        self.dt                 = dt #The time-step all filters are sampled at in seconds
        self.functionType       = [] #List of function types, for example [ ["gausLinSelf', a, b, c, timeLength, dt, timeDelay], ...]
        self.filterPositions    = [] #The list of indicies, where each index value points to a filterd neuron spike train stored in condInenFuncFilterManager.filteredArrays
        self.filterStepsDelayed = [] #The list of delays used for each of the filtered data.
        self.filterCoefficients = numpy.array([[]]) #The list of coefficients that go along with each filter.
        # This following array will change everytime we update data to be used to predict the states
        self.dataMatrix     = numpy.array([]) # The 2D matrix of training or input data data to the CIF, corresponds to the last data given.
        self.spikeArray     = numpy.array([]) # Boolean array
        self.spikeTimeArray = numpy.array([]) # Time-step values that correspond to spikes stored in self.spikeArray


    def doesFunctionExist(self, newFunctionType):
        """
        DESCRIPTION
        Indicates if the given "newFunctionType", which is a list of values that describes the function/filter, is
        alreayd specified for this particular CIF instantiation.

        :param newFunctionType: any of the acceptable lists for the self.functionType list, and are repeated below,
                                ['bias']
                                ['gausLinSelf', a, b, c, timeLength, dt, timeDelay]
                                ['gausLinSyn', otherNeuronGroup, otherNeuronId, a, b, c, timeLength, dt, timeDelay]
                                ['gausQuadSynSelf', otherNeuronGroup, otherNeuronId, a, b, d, timeLength, dt, timeDelaySelf, timeDelayOther]
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
        Find the index of the 'newFunctionType' in the self.functiionType variable.
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
        Add a new filter or functionType to the CIF object. NOTE, this method includes checking, thus before actually
        adding the new information, we make sure it is not allready included. If it already is, then the given info is
        NOT added as no duplicates are allowed.

        :param functionType:
        :param filterPosition:
        :param stepsDelayed:
        :return:
        """
        if not self.doesFunctionExist(functionType):
            self.functionType.append(functionType)
            self.filterPositions.append(filterPosition)
            self.filterStepsDelayed.append(stepsDelayed)
            if self.filterCoefficients.size>0:
                self.filterCoefficients = numpy.vstack((self.filterCoefficients,0.))
            else:
                self.filterCoefficients = numpy.array([[0.]])


    def appendSpikeArray(self, spikeTimes, startTime, endTime, dt):
        """
        DESCRIPTION
        For each neuron that has a CIF, we keep track of its spikes in memory, new spikes are appended here. In
        particular self.spikeArray and self.spikeTimeArray are appended. Data is only deleted if it overlaps currently
        stored values. We do not store spike-times explicity in this class, rather we have a boolean array (stored in
        self.spikeArray) and an array of equal size that indicates the time step (stored in self.spikeTimeArray)

        This method should be called when new simulation is added and filtered so that the data here is up to date.

        Note, methods will be provedied in this class to handle cleanup and deletion of memory in order to reduce memory
        overhead.

        :param spikeTimes:
        :param startTime:
        :param endTime:
        :return:
        """
        # Create boolean array in place of the spike times, remember we do not include the endTime time-step!
        timeArray = numpy.arange(start=startTime, stop=endTime, step=dt)
        boolArray = numpy.zeros(timeArray.size)
        indsArray = (spikeTimes-startTime)/dt
        # Make sure if there is a spike that has a spike-time of end-time, we make it slightly less so it rounds down
        indsArray[indsArray >= boolArray.size] -= 1
        boolArray[indsArray.astype(int)] = 1.
        # Append/update self.spikeArray and self.spikeTimeArray
        if self.spikeArray.size>0:
            # Find the next expected time-step that is thought to be seen occuring next in the input data. Depending on
            # on what this is will change how the appending process works
            expectedNextTime = self.spikeTimeArray[-1] + dt
            # Check if there is missing time between the stored and given data, thus we have to add filler time and
            # filler data. Then append the new information.
            if timeArray[0]>expectedNextTime:
                fillerTimes = numpy.arange(start=expectedNextTime, stop=startTime, step=dt)
                fillerZeros = numpy.zeros(fillerTimes.size)
                self.spikeArray = numpy.hstack((self.spikeArray, fillerZeros, boolArray))
                self.spikeTimeArray = numpy.hstack((self.spikeTimeArray, fillerTimes, timeArray))
            # Check if there is an overlap in time, if so we remove overlap and replace where needed.
            elif timeArray[0]<expectedNextTime:
                indsPre = numpy.where((self.spikeTimeArray<expectedNextTime))
                indsPost = numpy.where((self.spikeTimeArray>timeArray[-1]))
                if indsPre[0].size>0 and indsPost[0].size>0:
                    self.spikeArray = numpy.hstack((self.spikeArray[indsPre], boolArray, self.spikeArray[indsPost]))
                    self.spikeTimeArray = numpy.hstack((self.spikeTimeArray[indsPre], timeArray, self.spikeTimeArray[indsPost]))
                elif indsPre[0].size>0 and indsPost[0].size==0:
                    self.spikeArray = numpy.hstack((self.spikeArray[indsPre], boolArray))
                    self.spikeTimeArray = numpy.hstack((self.spikeTimeArray[indsPre], timeArray))
                elif indsPre[0].size==0 and indsPost[0].size>0:
                    self.spikeArray = numpy.hstack((boolArray, self.spikeArray[indsPost]))
                    self.spikeTimeArray = numpy.hstack((timeArray, self.spikeTimeArray[indsPost]))
                else:
                    self.spikeArray = boolArray
                    self.spikeTimeArray = timeArray
            # In this case, the input data directly follows (in time) the currently stored data, so this is a simple
            # appending process. This is the most desireable situation.
            else:
                self.spikeArray = numpy.hstack((self.spikeArray, boolArray))
                self.spikeTimeArray = numpy.hstack((self.spikeTimeArray, timeArray))
        # If the arrays are empty
        else:
            self.spikeArray = boolArray
            self.spikeTimeArray = timeArray


    def clearSpikeArray(self):
        """
        DESCRIPTION
        This is meant to clear out spike memory of the neuron in question. This should be typically run after
        probabilities are calculated and after CIF coefficients are updated.

        :return: N/A
        """
        self.spikeArray = numpy.array([])
        self.spikeTimeArray = numpy.array([])


    def createFeasibleFilteredDataSet(self, filteredArrays):
        """
        DESCRIPTION
        Here the filtered data is prepared to be used as data points for prediction of spike or no spike.

        We are parsing over 'filteredData' which is (or should be) litterally condInenFuncFilterManager.filteredArrays,
        which is the list of numpy.arrays([]), the filtered data spike data stored in the filter-manager.

        :param filteredArrays:
        :return: N/A
        """
        # We want to know how many time-steps that can be possibly be used for each of the filtered data arrays. In
        # order to do this we create an array 'maxSteps' that corresponds to each of the Cif storeded, and searches
        # through the filtered data to see how many time-steps we can keep.
        maxSteps = numpy.zeros(len(self.filterStepsDelayed))
        for ind, dlys in enumerate(self.filterStepsDelayed):
            # This can be less then zero in the case of a 'bias' term which will not have a filter to point to. In the
            # case of a 'bias' term, the value stored will be '-1', so we just set the value to some large value and go
            # with it
            if dlys>=0:
                maxSteps[ind] = filteredArrays[self.filterPositions[ind]].size - dlys
            else:
                maxSteps[ind] = 10000000
        maxStepsAllowed = int(numpy.amin(maxSteps)) # At most we take the minimum value.
        # Now given the maximum number of time-steps where we have a complete data set, we need to actually collect the
        # data and put into the correct form, in some N-D numpy.array([[]]) where each row is to be a sample of data,
        # and then fill the data matrix.
        if maxStepsAllowed > 0:
            self.dataMatrix = numpy.ones((maxStepsAllowed, len(self.filterPositions)))
            # Iterate over each column, which are the different data sources
            for i in range(self.dataMatrix.shape[1]):
                # If we have a bias term, we just keep the ones, otherwise we put the correct data in the arrays.
                if self.filterStepsDelayed[i]>0:
                    self.dataMatrix[:,i] = filteredArrays[self.filterPositions[i]][-(maxStepsAllowed+self.filterStepsDelayed[i]):-self.filterStepsDelayed[i]]
                elif self.filterStepsDelayed[i]==0:
                    self.dataMatrix[:,i] = filteredArrays[self.filterPositions[i]][-maxStepsAllowed:]
            # Now flip the matrix so that we have the lastest time along the bottom row, and the data from the earlier
            # times along the first row.
            self.dataMatrix = numpy.flipud(self.dataMatrix)
        else:
            # If there is no data, we just leave an empty matrix.
            self.dataMatrix = numpy.array([])


    def generateSpikeProbabilities(self):
        """
        DESCRIPTION
        Here the probabiliteis of a neuron spike is calculated. This calcualtion is based on data stored in
        self.dataMatrix which must have already be calculated. The data in self.dataMatrix is organized such that each
        row is data input for one time-step and each column represents a different component of the CIF.

        It is assumed that self.dataMatrix has been updated and given the most up to date data... and thus the data
        corresponds to the most up to date time values as well!

        :return:
        """
        # Numpy is weird, they use the method 'dot' to handle matrix multiplication.... ugg I do not like that
        if self.dataMatrix.shape[0]>0 and self.dataMatrix.shape[1]>0:
            ddd = 1
            glmOutput = numpy.dot(self.dataMatrix, self.filterCoefficients)
            glmOutput = 1./(1. + numpy.exp(-1.*glmOutput))
            dd = 2
        else:
            glmOutput = numpy.array([])
        return glmOutput


    def updateCifCoeffsSGD(self, learningRate):
        """
        DESCRIPTION
        Here, assuming both the self.dataMatrix and self.spikeArray are both up todate, we go ahead and run stochastic
        gradient descent to update the CIF coefficients.

        IMPORTANT note, it is assumed that the data stored in self.dataMatrix corresponds to the last values stored in
        self.spikeArray, essentially that the data was produced from the same time. This should be the case if the input
        data from the simulation is used to (1) update the filter-manager class obj, (2) then update the self.spikeArray
        and then (3) used to create the self.dataArray

        :param learningRate:
        :return:
        """
        if self.dataMatrix.shape[0]>0 and self.dataMatrix.shape[1]>0:
            # Generate the probability of a spike at each time-step according to the model currently stored.
            # predictionOutput = self.generateSpikeProbabilities()
            # Iterate over all the stored data points
            numbIts = self.dataMatrix.shape[0]
            for i in range(self.dataMatrix.shape[0]):
                predictionOutput = self.generateSpikeProbabilities()
                # Grab the last spikes (stored as 0s and 1s) in self.spikeArray based on how many we need (the number is
                # equal the size of the self.dataMatrix). Note that the last row in self.dataMatrix corresponds to the most
                # current time-step so we have to be careful about whic data point is grabbed.
                spikeVal = self.spikeArray[-(numbIts-i)]
                self.filterCoefficients = self.filterCoefficients + learningRate*(spikeVal - predictionOutput[i])*numpy.transpose(self.dataMatrix[i,:])


class condInenFuncFilterManager:
    """
    DESCRIPTION
    Class: manages different filter-funtions, filters spike trains, and stores results.

    General Guidelines on usage:
    (1) Adding new filters: In order to add a new filter-function for some neuron, one will typically first call
        self.addFilter, which will actually generate a filter (vector for convolution) and store the filter. Note, if
        the filter is already there it will not be duplicated. Note, if the filter is the same as other filters but with
        different delays, it is not duplicated either. The next typcial step is to call self.addNeuronAndFilter which
        will add the neuron and filter to and prepare the class to store filtered spike from that neuron. Again there
        is duplication protection. A number of other internal functions (self.getFilterIndex, self.getFilterTruncation,
        self.getFilterType, self.getNeuronAndFilterIndex) are used internally in the two formally mentioned methods.
    (2) Initial filter recorders: This class filters and stores the data. After all the filter-functions have been added
        then the user will typcially call self.initFilterRecorders(), which alots properly sized arrays in memroy for
        keeping parsed spike data in memory.
    (3) Filter spike data: After the memory for filtered spikes has been added, one can filter data by using the
        self.filterSpikeTrainData method.
    (4) Memory reduction, garbage collection: once the stored filtered data has been used for whatever purpose, for
        example prediction spikes or updating CIF coefficients, then it can be removed for memory consolidation. The
        most appropriate methods to call are self.garbageCollectFilteredArrays() and self.garbageCollectSpikeTrainMemory
        which eliminate most but not all data. All the appropriate data is stored in memory such that all filters that
        look at delayed time-steps will be operable in the next time-step.
    """

    def __init__(self):
        """
        DESCRIPTION
        Conditional Intensity Filter Manager, this class manages all of the filters used to construct the CIFs for each
        given neuron. In this class there are 3 varaibles.

        self.filterIdentities: list of list, [ [filter 1 info], [filter 2 info], ... ], where in each sublist the
                               following variants are acceptable,
                               ['gaussian', a, b, c, timeLength, dt, timeDelay, truncation, filterIndex]
        self.filters:          list of 1D numpy.array, filters where each filter is something we convolve the spike-
                               trains that have been recorded.
        self.filtersDt:        list of time-step values for each filter, this list should be as long as the
                               self.filters list.
        self.filtersBackUpTime list of the amount of time we want to inlcude backwards in order to make the filters
                               accurate. This is because we do not want to neglect spikes that were formallay recorded.
        self.neuronToFilter    2D numpy.array, an example row being [ [neuron-group, neuron Id, filterIndex], ... ],

        self.spikeTrainMemory  a list of 1D numpy.array, where in each array we have all the spike times, [ [t1, t2, ...
                               tn], ... ]
        self.filteredArrays    a list of 1D numpy.array, where in eacy array we have all the
        """
        #The size of the following list is based on the number of unique fiters and nuron combos
        self.filterIdentities     = [] #Association between filter-identity/filter-type to an actual filter in self.filters

        #The size of the following lists is eqaul to the number of unique convolution filters needed
        self.filters              = [] #The actual convoltuion filters which are 1D numpy.arrays
        self.filtersDt            = [] #The time-step associate to filter in self.filters
        self.filtersBackUpTime    = [] #The amount of time before we need to include backwards for each filter in self.filters

        #The size of the following are equal to the number of unique neuron-convoFilter pairs
        self.neuronToFilter            = numpy.array([[]]) # The association between filters and neurons, numpy.array([ [neuron-group, neuron Id, filterIndex], ...]
        self.numbStepsToStore          = numpy.array([[]]) # The number of time-steps we need to keep for the next time we run calculations, essentially this is the maximum delay we need
        self.actualDelay               = numpy.array([[]]) # Realy just kept as a measure to debug, this tells us how much of a delay each filter has.

        self.filteredArrays            = [] # The filtered spike-train output at each time-step resultant from convolution operation
        self.filteredArraysTimeValues  = [] # The time-step corresponding to each filtered output time-step, should be the same size as the filterArray
        self.spikeTrainMemory          = [] # The old spike times we want to keep for each neuron


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


    def getFilterType(self, filterIndex):
        """
        DESCRIPTION
        This function is primarily for debugging. A 'filterIndex' is given, which correspond to the index of
        (1) self.filters, (2) self.filtersDt, and (3) self.fitlersBackUpTime. Given this number we want to see which
        function-type describes the filter stored in self.filters[fitlerIndex]

        :param filterIndex: integer, the index of a particular filter in self.filters, (self.fitlers is a list of 1D
                            numpy.arrays() that are convolved with a spike train).

        :return: list of filter information, in the case the filter is of 'gaussian' type, then we return a list with
                 ['gaussian', a, b, c, timeLength, dt]
        """
        toReturn = ["NO FUNCTION STORED"]
        for filtTypes in self.filterIdentities:
            if filterIndex == filtTypes[-1]:
                if filtTypes[0] == 'gaussian':
                    toReturn = filtTypes[0:6]
                    break
        return toReturn


    def getNeuronAndFilterIndex(self,neuronGroupId, neuronId, filterIndex):
        inds = numpy.where((self.neuronToFilter == (neuronGroupId, neuronId, filterIndex)).all(axis=1))
        if inds[0].size > 0:
            dummyOut = inds[0][0]
        else:
            dummyOut = -1
        return dummyOut


    def addFilter(self, filterInfo):
        """
        DESCRIPTION
        Adds a new filter (a filter in this case means something specifically that is convolved with spike trains) if
        the filter described by the input 'filterInfo' does not alreayd exist. Imporanatly, this method includes internal
        checks to make sure no duplicates are added. Thus, if a filter is specified in the input that is already stored,
        a dumplicate will NOT be added.

        An example of 'filterInfo is ['gaussian', a, b, c, timeLength, dt, timeDelay]. Now remember, there maybe more
        self.filterIdentities than actual self.fitlers. This is because many of those filters described in the
        self.filterIdentities are actually produce the same output just delayed by different time-steps!

        :param filterInfo: a list of filter information, where the following varients are acceptable,
                           ['gaussian', a, b, c, timeLength, dt, timeDelay]

        :return: [truncation, filterIndex, timeStepsDelayed, implicitDelay]
        """
        if filterInfo[0]=='gaussian':
            # First determine the number of time-steps we need to look backwards in memory of the filtered spike-trains.
            # That is some filters are time-delayed gaussians. Each of the fitlers created implecitly have a delay
            # immbedded. However if the the the time-dealy is longer then the filter's 'timeLength', then we actually
            # have to store a dealy value!.
            timeDelay  = filterInfo[6]
            timeLength = filterInfo[4]
            truncation = min(timeLength, timeDelay)
            timeStepsDelayed = int((timeDelay - truncation)/filterInfo[5])
            if self.getFilterIndex(filterInfo[0:6], timeDelay=filterInfo[6]) != -1:
            #if self.getFilterIndex(filterInfo[0:6], timeDelay=None, truncation=truncation):
                # This filter already exists, which means we only have to grab the relevant output details which in this
                # case is the size of the filter's truncation, and the filter-index.
                dummyList = [self.getFilterTruncation(filterInfo[0:6],timeDelay=filterInfo[6]), self.getFilterIndex(filterInfo[0:6],timeDelay=filterInfo[6]), timeStepsDelayed]
            else:
                # Filter of that type does not exist, Remember, though some filters are different,
                # they essentially are the same filter convolved with a spike train, but with differnet delay values.
                if self.getFilterIndex(filterInfo[0:6], truncation=truncation) != -1:
                    # Filter exists with different time-delay (but really same filter/convolution operation). Record the
                    # filterInfo, but do not create a new filter for convolution.
                    newFilterIndex = self.getFilterIndex(filterInfo[0:6], truncation=truncation)
                    filterInfo.append(truncation) #Append to the filter information the truncation time
                    filterInfo.append(newFilterIndex) #Append to the filter the index of the filter in self.filters
                    self.filterIdentities.append(filterInfo)
                else:
                    # In this case we have to create a new filter because it does not already exist. Firstly we want
                    # to determine the filter-backup-time. This is the amount of time over some backwards horizon
                    # wherein we need to store neuron spikes in memory so the filter operations are accurate.
                    filBackUpTime = truncation - (-1.*timeLength)
                    self.filters.append(truncatedGaussianFilterGenerator(filterInfo[1], filterInfo[2], filterInfo[3], filterInfo[4], filterInfo[5], truncation))
                    self.filtersDt.append(filterInfo[-2])
                    self.filtersBackUpTime.append(filBackUpTime)
                    newFilterIndex = len(self.filters)-1
                    filterInfo.append(truncation)
                    filterInfo.append(newFilterIndex)
                    self.filterIdentities.append(filterInfo)
                dummyList = [truncation, newFilterIndex, timeStepsDelayed]
        elif filterInfo[0]=='bias':
            # We do not add 'bias' filters as with this type of CIF variable, there is no spike train filtering at all
            dummyList = [-1, -1, -1]
        return dummyList


    def addNeuronAndFilter(self, neuronGroupId, neuronId, filterIndex, timeStepsDelayed, timeDelay):
        """
        DESCRIPTION
        This class contians/stores filters (which are things we convolve spike trains with). It also contains variables
        to store the results of the convoltuion with spike-trains. The question that remains, what filters need to be
        applied to which spike trains? This problem is solved with self.neuronsToFilter which is a 2D numpy.array of
        size NX3, where each row has the format [neuron-group, neuron id, filterIndex]. The filter index in this case
        is the index of the filter stored here in self.filters.

        Here a new row is added to self.neuronsToFilter only if, that row does NOT already exist. That is, internatlly
        this method will check to make sure this is not a duplicate entry. If it is, then a new row will not be added.

        :param neuronGroupId: integer, the neurons's group id
        :param neuronId: interger, the neurons index in its group
        :param filterIndex: the index in self.neuron
        :param timeStepsDelayed: how many steps backwards in filterd data we have to look, some filters have implicit
                                 delays and will not actually have any steps to look backwards.
        :return: the index of the self.filteredArrays we will access to get the data needed for CIFs
        """
        if self.neuronToFilter.size == 0: # It the original array is empty, then we need to make the first row
            self.neuronToFilter = numpy.array([[neuronGroupId, neuronId, filterIndex]])
            self.numbStepsToStore = numpy.array([[timeStepsDelayed]])
            self.actualDelay = numpy.array([[timeDelay]])
            dummyOut = 0
        else:
            inds = numpy.where((self.neuronToFilter==(neuronGroupId, neuronId, filterIndex)).all(axis=1))
            if inds[0].size > 0:
                dummyOut = inds[0][0] #This is the position of the filtered values we are going to want to sample
                # Check if we need to update the number of time-steps we need to store.
                if timeStepsDelayed > self.numbStepsToStore[inds[0][0]]:
                    self.numbStepsToStore[inds[0][0]] = timeStepsDelayed
            else:
                # The filter and neuron combo does not exist lets add it
                self.neuronToFilter = numpy.vstack((self.neuronToFilter, numpy.array([neuronGroupId, neuronId, filterIndex])))
                self.numbStepsToStore = numpy.vstack((self.numbStepsToStore, numpy.array([timeStepsDelayed])))
                self.actualDelay = numpy.vstack((self.actualDelay, numpy.array([timeDelay])))
                dummyOut = self.neuronToFilter.shape[0]-1 #Must subtract the 1 because we do zero indexing
        return dummyOut


    def initFilterRecorders(self):
        """
        DESCRIPTION
        Once all the filters needed have been given via the class method self.addFilter(), and all the neurons that need
        their spike trains filtered have been indicates via the class method self.addNeuronAndFilter(), then this
        function should be called.

        This function initializes the class variables that store all the spikes and filtered spiek trains.

        It is only allowed to call this operator once. Every time it gets called if will overwrite all formally stored
        data.

        NOTE: some of the memory will be duplicates (at the moment) where we are repeatedly storing the spikes from the
        same neuron over and over again.

        :return:
        """
        #Simply initializing here, this WILL overwrite data if it already exists
        self.filteredArrays = [numpy.array([]) for i in range(self.neuronToFilter.shape[0])]
        self.filteredArraysTimeValues = [numpy.array([]) for i in range(self.neuronToFilter.shape[0])]
        self.spikeTrainMemory = [numpy.array([]) for i in range(self.neuronToFilter.shape[0])]


    def filterSpikeTrainData(self, listOfSpikeObjs, listOfNeuronGroups, startTime, endTime):
        """
        DESCRIPTION
        Filter though given spike-trains. The given spike trains must be of the form brian2.SpikeMonitor class objects.
        These brian2.SpikeMonitor objects are given in the list 'listOfSpikeObjs'. We know which neuron-group each of
        the brian2.SpikeMonitor objects are associated with by the also given input 'listOfNeuronGroups'. These two
        inputs must be lists of the same length as they correspond to eachother.

        It is important to make sure that 'startTime' is far enough backwards that we have an overlap in time so that
        the convolution filters are accurately calculated.

        IMPORTANT to note, we assume data is given each time for ALL neurons that are being tracked in
        self.neuronsToFilter! This is very important. If there is no spike data for one of those neurons, it is assumed
        that said neuron produced zero spikes from startTime to endTime.

        IMPORTANT to note, the last time-step at time=endTime is NOT included or stored.

        IMPORTANT to note, the brain2.SpikeMonitor objects must index their neurons, in the same way that the neurons
        have been index in all the filter-manager variables. For example ,if a neuron belongs to neuron-group 2, with
        neuron-id of 5, then we will expect the brain2.SpikeMonitor object for neuron-group2 to have the spikes stored
        in said object ot be associated with index 5.

        :param listOfSpikeObjs:
        :param listOfNeuronGroups:
        :param startTime:
        :param endTime:
        :return:
        """
        # Iterate through neurons that need to be filtered, as noted in this class's self.neuronToFilter numpy.arrray.
        # That is we iterate over all the neurons that need to be filtered and look to see if they are inlcuded in the
        # given data. Remember self.neuronToFilter is a 2D numpy.array where each row is
        # [neuron-group, neuron Id, filterIndex]
        for ind in range(self.neuronToFilter.shape[0]):
            nGrp  = self.neuronToFilter[ind, 0] # Integer, neuron's group Id
            nId   = self.neuronToFilter[ind, 1] # Integer, neuron's inter-group Id
            filt  = self.neuronToFilter[ind, 2] # Integer, filter index (references self.filters, self.filtersDt, self.filtersBackUpTime, etc.)
            dStps = self.numbStepsToStore[ind]  # Integer, the number of time-steps we keep in memroy because of the delays with some filters.
            # ONLY the neurons are going to be updated if they have input data to parse through
            if nGrp in listOfNeuronGroups: #Check if the input data has the neuron group represented or not
                spkObjInd = listOfNeuronGroups.index(nGrp) #Indicates which of the brian2.SpikeMonitor objects to grab from 'listOfSpikeObjs' list
                spkObj = listOfSpikeObjs[spkObjInd] #Grab the correct brian2.SpikeMonitor object
                # 'spkObj' is a brian2.SpikeMonitor object. It may have spike times that preccede 'startTime' (here
                # given as input. As such we need to ONLY look at consider those times >=startTime, and <=endTime, and
                # lastly that we are looking at the neuron in question we care about.
                timeCorrectedIndsObj = numpy.where((spkObj.t>=startTime) & (spkObj.t<=endTime) & (spkObj.i==nId)) #Have to use the '&' bitwise operation here .. not really sure why yet
                timeCorrectedIndsMem = numpy.where( (self.spikeTrainMemory[ind]>=startTime-self.filtersBackUpTime[filt]) & (self.spikeTrainMemory[ind]<=startTime) ) #Have to use the '&' bitwise operation here
                if timeCorrectedIndsObj[0].size>0 and timeCorrectedIndsMem[0].size>0:
                    spikeTimes = numpy.hstack((spkObj.t[timeCorrectedIndsObj[0]], self.spikeTrainMemory[ind][timeCorrectedIndsMem[0]]))
                elif timeCorrectedIndsObj[0].size>0 and timeCorrectedIndsMem[0].size==0:
                    spikeTimes = spkObj.t[timeCorrectedIndsObj[0]]
                elif timeCorrectedIndsObj[0].size==0 and timeCorrectedIndsMem[0].size>0:
                    spikeTimes = self.spikeTrainMemory[ind][timeCorrectedIndsMem[0]]
                else:
                    # No spikes given as input, and no spikes in memory, we assume no spikes over the time-period in
                    # question and will update the rest as such.
                    spikeTimes = numpy.array([])
                # Calculate filtered spike train data, and generate the associted time-array, both of which will be
                # 1D numpy.array objects. Then update the the appropriate varibles in this class. NOTE, the
                # convolution method will NOT include the last time-step!
                filteredSpikeTrainList = convolveSpikeTimesWthBackUp(self.filters[filt], spikeTimes, startTime, endTime, self.filtersDt[filt], backupTime=self.filtersBackUpTime[filt])
                # Call functions to append filtered data, and spike data to memory.
                self.appendFilteredArrays(ind, filteredSpikeTrainList[0], filteredSpikeTrainList[1])
                self.appendSpikeTrainMemory(ind, spikeTimes)


    def appendFilteredArrays(self, ind, filteredArray, filteredArrayTimes):
        """
        DESCRIPTION
        A method, called generally when given new spike data that has been filtered, to append self.filteredArrays
        and self.filteredArraysTimeValues. Here NO DATA is deleted, only data is added. Data should only be deleted
        once it is used to calculate CIFs or update CIFs. This step comes prior to those steps.

        IMPORTANT to note, data is deleted if there is an overlap in time-arrays

        IMPORTANT to note, if there is a gap in time between the stored memory and the given data, that time is filled
        in with blank data to make there is no gap in time.

        :param ind: integer, the index of the current neuron / filter combination stored in self.neuronToFilter where
        :param filteredArray: numpy.array([])
        :param filteredArrayTimes: numpy.array([])
        :return: N/A
        """
        if filteredArray.size>0: # Check if the input data is valid
            #Find the time-step!
            dt = self.filtersDt[self.neuronToFilter[ind,2]]
            # Check if the stored array previously has values.
            if self.filteredArraysTimeValues[ind].size>0:
                # Find the next expected time-step that is thought to be seen occuring next in the input data. Depending
                # on what this is will change how the appending process works.
                expectedNextTime = self.filteredArraysTimeValues[ind][-1] + dt
                # Check if there is missing time between the stored and given data, thus we have to add filler time and
                # filler filtered data. Then append the correct self.filtereArrays and self.filterArraysTimeValues
                if filteredArrayTimes[0]>expectedNextTime:
                    fillerTimes = numpy.arange(start=expectedNextTime, stop=self.filteredArraysTimeValues[ind][0], step=dt)
                    fillerZeros = numpy.zeros((1,fillerTimes.size))
                    self.filteredArrays[ind] = numpy.hstack((self.filteredArrays[ind], fillerZeros, filteredArray)) #Filtered spike train
                    self.filteredArraysTimeValues[ind] = numpy.hstack((self.filteredArraysTimeValues[ind], fillerTimes, filteredArrayTimes)) #Associated time arrays
                # Check if there is an overlap in time, between stored data and that of the newly given data. This is
                # not desireable, best is that given data is a continuation of some simulation. However, if this is the
                # case, we clear out part of the memory in question, and replace it. It is UP TO THE USER to make sure
                # that this operation makes sense.
                elif filteredArrayTimes[0]<expectedNextTime:
                    indsPre = numpy.where((self.filteredArraysTimeValues[ind]<expectedNextTime))
                    indsPost = numpy.where((self.filteredArraysTimeValues[ind]>filteredArrayTimes[-1]))
                    if indsPre[0].size>0 and indsPost[0].size>0:
                        self.filteredArrays[ind] = numpy.hstack((self.filteredArrays[ind][indsPre], filteredArray, self.filteredArrays[ind][indsPost]))
                        self.filteredArraysTimeValues[ind] = numpy.hstack((self.filteredArraysTimeValues[ind][indsPre], filteredArrayTimes, self.filteredArraysTimeValues[ind][indsPost]))
                    elif indsPre[0].size>0 and indsPost[0].size==0:
                        self.filteredArrays[ind] = numpy.hstack((self.filteredArrays[ind][indsPre], filteredArray))
                        self.filteredArraysTimeValues[ind] = numpy.hstack((self.filteredArraysTimeValues[ind][indsPre], filteredArrayTimes))
                    elif indsPre[0].size==0 and indsPost[0].size>0:
                        self.filteredArrays[ind] = numpy.hstack((filteredArray, self.filteredArrays[ind][indsPost]))
                        self.filteredArraysTimeValues[ind] = numpy.hstack((filteredArrayTimes,self.filteredArraysTimeValues[ind][indsPost]))
                    else:
                        self.filteredArrays[ind] = filteredArray
                        self.filteredArraysTimeValues[ind] = filteredArrayTimes
                # In this case, the input data falls directly following the currently stored data, so this is a simple
                # appending process. This is the most desireable situation where new data given is immediatly follows
                # the formally given data.
                else:
                    self.filteredArrays[ind] = numpy.hstack((self.filteredArrays[ind], filteredArray))  # Filtered spike train
                    self.filteredArraysTimeValues[ind] = numpy.hstack((self.filteredArraysTimeValues[ind], filteredArrayTimes))  # Associated time arrays
            # If self.filteredArrays is empty (nothing yet stored) there is no special operation, just saving to the
            # correct variable in memory.
            else:
                self.filteredArrays[ind] = filteredArray
                self.filteredArraysTimeValues[ind] = filteredArrayTimes


    def garbabgeCollectFilteredArrays(self):
        '''
        DESCRIPTION
        If this method is called, it is assumed that all the CIF prediction and updating for some epoch has concluded.
        This method will go through and shrink the numpy.arrays in self.filteredArrays and in
        self.filteredArraysTimeValues to the smallest sizes possible as denoted by the values stored in
        self.numbStepsToStore

        :return: N/A
        '''
        # Itearte over all the neuron / filter combinations saved in self.neuronToFilter
        for ind in range(self.numbStepsToStore.size):
            dStps = self.numbStepsToStore[ind] # The number of time-steps we want to save in the filtered arrays.
            if dStps>0:
                self.filteredArrays[ind] = self.filteredArrays[ind][-dStps:]
                self.filteredArraysTimeValues[ind] = self.filteredArraysTimeValues[ind][-dStps:]
            else:
                self.filteredArrays[ind] = numpy.array([])
                self.filteredArraysTimeValues[ind] = numpy.array([])


    def clearFilteredArrays(self, ind, startTime=0., stopTime=-1):
        """
        DESCRIPTION
        A method to clear diltered arrays, for a particular neuron / filter indicated by the postition in
        self.neuronToFilter by the input 'ind', over some time-range given by starTime - infinity, or startTime to
        stopTime

        :param ind:
        :param startTime:
        :param stopTime:
        :return:
        """
        if stopTime < startTime:
            timeCorrectedInds = numpy.where((self.filteredArraysTimeValues[ind] >= startTime))
            if timeCorrectedInds[0].size > 0:
                self.filteredArraysTimeValues[ind] = numpy.delete(self.filteredArraysTimeValues[ind],
                                                                  timeCorrectedInds[0])
                self.filteredArrays[ind] = numpy.delete(self.filteredArrays[ind], timeCorrectedInds[0])
        else:
            timeCorrectedInds = numpy.where(
                (self.filteredArraysTimeValues[ind] >= startTime) & (self.filteredArraysTimeValues[ind] <= stopTime))
            if timeCorrectedInds[0].size > 0:
                self.filteredArraysTimeValues[ind] = numpy.delete(self.filteredArraysTimeValues[ind],
                                                                  timeCorrectedInds[0])
                self.filteredArrays[ind] = numpy.delete(self.filteredArrays[ind], timeCorrectedInds[0])


    def appendSpikeTrainMemory(self, ind, spikeTimes):
        """
        DESCRIPTION
        A method, called generally when given new spike data that has been filtered, to append self.spikeTrainMemory
        Here data is only deleted, if the new formerly stored spikes overlap with the new range of time that the new
        set of spikes given in 'spikeTimes'. Otherwise, the spikes are simply appended.

        :param ind:
        :param spikeTimes:
        :return:
        """
        if spikeTimes.size>0: # Check input
            if self.spikeTrainMemory[ind].size>0:
                # Remove all those times that overlap with newly given spike times if such spikes exist.
                firstNewSpikeTime = numpy.amin(spikeTimes)
                lastNewSpikeTime = numpy.amax(spikeTimes)
                self.spikeTrainMemory[ind] = self.spikeTrainMemory[ind][ (self.spikeTrainMemory[ind]<firstNewSpikeTime) & (self.spikeTrainMemory[ind]>lastNewSpikeTime) ]
                if self.spikeTrainMemory[ind].size>0:
                    # Add the new spikes to our spike train
                    self.spikeTrainMemory[ind] = numpy.hstack((self.spikeTrainMemory[ind], spikeTimes))
                else:
                    self.spikeTrainMemory[ind] = spikeTimes
            else:
                self.spikeTrainMemory[ind] = spikeTimes
            # Sort the spike times!
            self.spikeTrainMemory[ind] = numpy.sort(self.spikeTrainMemory[ind])


    def garbageCollectSpikeTrainMemory(self, endTime):
        '''
        DESCRIPTION
        If this method is called, it is assumed that all the CIF prediction and updating for some epoch has concluded.
        This method will go through and shrink the numpy.arrays in self.spikeTrainMemory to the smallest sizes possible
        as denoted by ... the value stored in self.filtersBackUpTime

        :return: N/A
        '''
        for ind, bckTime in enumerate(self.filtersBackUpTime):
            self.spikeTrainMemory[ind] = self.spikeTrainMemory[ind][ (self.spikeTrainMemory[ind]<endTime) & (self.spikeTrainMemory[ind]>endTime-bckTime) ]


    def clearSpikeTrainMemory(self, ind, startTime=0., stopTime=-1):
        """
        DESCRIPTION
        A method to delete spike times, for a particular neuron / filter indicated by the postition in
        self.neuronToFilter by the input 'ind', over some time-range given by starTime - infinity, or startTime to
        stopTime

        :param ind: the neuron index, essentially the index that lines up with self.neuronToFilter
        :param startTime: the time beginning from which we want to delete stored spikes
        :param stopTime:  the time to which we want to delete spikes, a value < 0 indicates to delete all spikes greater
                          in time then 'startTime'
        :return:
        """
        if stopTime < startTime:
            timeCorrectedInds = numpy.where((self.spikeTrainMemory[ind] >= startTime))
            if timeCorrectedInds[0].size > 0:
                self.spikeTrainMemory[ind] = numpy.delete(self.spikeTrainMemory[ind], timeCorrectedInds[0])
        else:
            timeCorrectedInds = numpy.where(
                (self.spikeTrainMemory[ind] >= startTime) & (self.spikeTrainMemory[ind] <= stopTime))
            if timeCorrectedInds[0].size > 0:
                self.spikeTrainMemory[ind] = numpy.delete(self.spikeTrainMemory[ind], timeCorrectedInds[0])



class condItenFuncManager:
    """
    DESCRIPTION
    Class: manages all the conditional intensity functions for all the given neurons.

    KNOWN PROBLEMS:
    (1) july-7-2016: I think the way I give indicate start and stop time does not make sense... needs to be checked. In
                     particular, when I give tStart and tEnd, currently I include tEnd. However this DOES NOT make
                     sense, as all teh spike would have occured before tEnd, and thus would be included ... right??
                     or do I not include the start time???
                     The answer is going to be I do not include the end-time, that is if I start simulation at t=0, and
                     then from for one dt=.001ms, then I will only have data from 1-time-step which is annotated as
                     time-step 0,
                     ... this might be fixed, need to run and test to be sure
    """


    def __init__(self, dt):
        """
        DESCRIPTION
        This is a class that manages all of the conditional intensity functions for all the neurons.
        """
        # An instantiattion of the filter-manager class
        self.filterManager = condInenFuncFilterManager()
        # Where we keep all our individual conditional intensity function objects for each individual neurons
        self.condIntenFunctions = []
        # The time-step we expect all the filters to be sampled at!!  I NEED TO REWRITE A BUNCH OF STUFF SO THAT WHEN
        # NEW FILTES ARE GIVEN, dt IS NO LOGER CONSIDERED OR NEEDED, RATHER WE SIMPLY REFERENCE THE GLOBAL VALUES HERE!
        self.dt = dt


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


    def createNeuronCifObj(self, neuronGroup, neuronId, dt):
        self.condIntenFunctions.append(conItenFunc(neuronGroup, neuronId, dt))
        return self.condIntenFunctions[-1] #Return the last created cifObj, which is the one we last created.


    def addNewFilter(self, newFilter):
        '''
        DESCRIPTION
        Adds a newfilter to a perscribe set of neurons. More specifically, the input (see input description below)
        describes a type of filter/function for a CIF, for a neuron or neurons. The class here then alters or addusts
        its list of CIF objects (self.condIntenFunctions) and its filter-manager object (self.filterManager). Quickly
        we will outlien the process,
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
        :return: Nothing retuned
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
                                    # Ammend the 'cifObj' so we know which filter to point to.
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
                        outVals0 = self.filterManager.addFilter( ['gaussian', newFilter[3], newFilter[4], newFilter[4],
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
        This method is (1) simply pushes the input data the filter-manager class (self.filterManger), which will
        automatically filter the data, and store it in the correct form, and (2) indicate to the CIF objects that new
        data is available to to used to predict spikes and used to train and update CIR coefficients.

        :param listOfSpikeObjs:
        :param listOfNeuronGroups:
        :param startTime:
        :param endTime:
        :return:
        """
        # FIRST, give the data to the spike train-manager to let it parse through the data, and filter what it needs to.
        self.filterManager.filterSpikeTrainData(listOfSpikeObjs, listOfNeuronGroups, startTime, endTime)
        # SECONDLY, now the spikes of the neurons need to be reflected in each of the CIF objects. This is done by by
        # iterating over the conditional intensity functions objects stored in the the list, self.condIntenFunctions
        for cifObj in self.condIntenFunctions:
            # Determine if and where our neurons are represented.
            if cifObj.neuronGroup in listOfNeuronGroups:
                spkObjInd = listOfNeuronGroups.index(cifObj.neuronGroup)
                spkObj = listOfSpikeObjs[spkObjInd]
                # 'spkObj' is a brian2.SpikeMonitor object. It may have spike tiems that preceede 'startTime'. As such
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


    def prepareParsedData(self):
        """
        DESCRIPTION
        This method gives the fitlered data stored in the filter-manager and gives it to the CIF objects such that they
        can ... predict spikes
        :return:
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

        :param iterations: how many times do we want to iterate over the given data
        :param method:
        :param learningRate:
        :return:
        """
        for cifObj in self.condIntenFunctions:
            cifObj.updateCifCoeffsSGD(learningRate)


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
                print("... Neuron",i.neuronId," of Neuron-Group",i.neuronGroup)
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
            plt.plot(times,vals)
        plt.show()


    def diagnosticCheckExamineConstrutedDataSet(self):
        print("\n")
        print("#######################################################################################################")
        print("#                                                                                                     #")
        print("#                                    Examining Constructred Data                                      #")
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
                #Note the data stored in the 'dataMatrix' numpy.array is such that the earlest time is given first, and
                #the lastest time is last, we need to flip the data so it shows up correctly
                plt.plot(cifObj.spikeTimeArray[-cifObj.dataMatrix.shape[0]:], numpy.flipud(cifObj.dataMatrix[:,columnInd]))
            plt.plot(cifObj.spikeTimeArray[-cifObj.dataMatrix.shape[0]:], cifObj.spikeArray, '--', linewidth=3.0)
        plt.show()







########################################################################################################################
# Internal Unit Testing
########################################################################################################################


if __name__ == "__main__":


    #################Testing the filter manager bits
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


    #################Testing the addition, and filtering of data
    # Create a new conditional-intensity-function manager object
    myFilterMan = condItenFuncManager(.001)
    # Create neurons and neurons groups
    neuronGroup0Ids = numpy.array([0])
    neuronGroup1Ids = numpy.array([0,1])
    #synapses0       = numpy.array([[0,1,0,1],[0,1,1,1]])
    synapses0       = numpy.array([[0,1,1,1]])
    gausFiltDelays  = numpy.array([0., .005, .01, .015, .02])
    # Create new conditional intensity functison for neurons
    inputFunctionType = ['gausLinSyn', 1, neuronGroup1Ids, synapses0, 1., .0, .002, .1, .001, gausFiltDelays] #['gausLinSyn', neuronGroup, neuronIds, synStruc, a, b, c, timeLength, dt, timeDelays]
    myFilterMan.addNewFilter(inputFunctionType)
    # Create new conditional intensity functison for neurons
    inputFunctionType = ['bias', 1, neuronGroup1Ids]
    myFilterMan.addNewFilter(inputFunctionType)
    # Check the filters
    myFilterMan.diagnosticCheckNeuronAndFilterCombs()
    #myFilterMan.diagnosticCheckExamineStoredFilters()
    #Now I create and add some data
    class someSpikeObj:
        def __init__(self, inds, times):
            self.t = times
            self.i = inds
    neuSpks = numpy.array([0,0,0,0,0,0,0])
    neuTmes = numpy.array([.001, .0081, .011, .025, .050, .0493, .08])
    neuTmes = neuTmes + .0001
    grp1SpikeMon = someSpikeObj(neuSpks,neuTmes)
    neuSpks = numpy.array([0, 1, 0, 1])
    neuTmes = numpy.array([.011, .025, .06, .075])
    neuTmes = neuTmes + .0001
    grp2SpikeMon = someSpikeObj(neuSpks,neuTmes)
    # Initialize the filter class
    myFilterMan.filterManager.initFilterRecorders()
    # Now actually give the data to be presented
    myFilterMan.parseNewSpikeTrainData([grp1SpikeMon,grp2SpikeMon], [0,1], startTime=0., endTime=.125)
    # Now some debugging
    #myFilterMan.diagnosticCheckExamineFilteredSpikeData()
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
    #myFilterMan.diagnosticCheckExamineConstrutedDataSet()
    # Now we want to try and train using gradient descent!
    def testSgdAndPrections(myFilterMan, trainIterations, cifObjNumb):
        # Grab the cif-object
        cifObj = myFilterMan.condIntenFunctions[cifObjNumb]
        # Begin to plot things by plotting out the original data!
        plt.figure(1)
        plt.plot(cifObj.spikeTimeArray[-cifObj.dataMatrix.shape[0]:], cifObj.spikeArray[-cifObj.dataMatrix.shape[0]:], '--', linewidth=3.0)
        # Generate initial predictions, and plot them
        probs = cifObj.generateSpikeProbabilities()
        plt.plot(cifObj.spikeTimeArray[-probs.size:],probs)
        #Iteritively train using SGD method
        for i in range(trainIterations):
            cifObj.updateCifCoeffsSGD(.01)
            probs = cifObj.generateSpikeProbabilities()
            print(probs.shape)
            plt.plot(cifObj.spikeTimeArray[-probs.size:], probs)
        plt.show()
    # Run the test method we just created
    testSgdAndPrections(myFilterMan, 5, 1)








