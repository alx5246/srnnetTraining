# A. Lonsberry
# September 2016 - Taken from sandbox (condIntenFunc_0_3.py)
#
# DESCRIPTION
#   The class(es) below are designed to represent and hold all the relevant data for a single neuron in a network.
#   That is this class will handle calculating and finding CIF output and learning for a single neuron.

import numpy
import matplotlib.pyplot as plt
import numpy.matlib
import time

# Also import comiled Code
#import updateCifCoeffsSGD
from cythonBuilt import updateCifCoeffsSGD #When compiled code is saved within a seperate directory

class condItenFunc:
    """
    DESCRIPTION
    A class that contains all the conditional-intensity function information for a single neuron. More specifically,
    this stores a list of off the individual functions that comprise the CIF (linear combinations constitute the entire
    CIF) and stores an array indicating which filtered data in some condInenFuncFilterManager (see class below) to get
    data from.

    GENERAL USE GUIDELINES:
    (1) Add functions to the CIF object: typically after an object of this class is created, you can add the functions
        that comprise the GLM suing the self.addFunction() method. A CIF as it is here described is a GLM that is passed
        through a logistic function. Every time one uses self.addFunction() method you are adding another linear
        component to the GLM

    (2) Add spikes and create data: once a simulation has been run and data is available there are two methods that need
        to be called. First is the self.appendSpikeArray() which takes in the new spikes and adds them to the memory
        here. Second is the self.createFeasibleFitleredDataSet() which takes in data from  a filter-manager object and
        creates a matrix of data that can be used to generate probabilities.

    (3) Find probabilities of spikes using data crated by running self.createFeasibleFilteredDataSet() to create data
        and thrn running self.generateSpikeProbabilities() to actually crate probabiliteis.

    (4) Update the GLM coefficients by using the stored data and created data generated by
        self.createFeasibleFilteredDataSet()

    THINGS TO ADD
    (1) JULY 2016 : Garbage collection. Like in my other methods I need a method to garbage collect on saved spikes
        that are in memory.
    """


    def __init__(self, neuronGroup, neuronId, dt):
        """
        self.neuronGroup: integer, the group in which the neuron is from
        self.neuronId: integer, the neuron index or ID associated to the neuron within some group
        self.functionType: list of lists for the following forms,
            ['bias']
            ['gausLinSelf', a, b, c, timeLength, dt, timeDelay]
            ['gausLinSyn', otherNeuronGroup, otherNeuronId, a, b, c, timeLength, dt, timeDelay]
            ['gausQuadSynSelf', otherNeuronGroup, otherNeuronId, a, b, d, timeLength, dt, timeDelaySelf, timeDelayOther]
        self.filterPoistions: list of indices, where each index stored points to one of the indices of filtered spike
            data that is stored in some condInenFucFilterManager object (see condInenFuncFilterManager.filteredArrays),
            The condInenFuncFilterManager.fitleredArrays is a list of 1D numpy.array objects. Each of the these objects
            is the result of convolving some function (indicated appropriately here self.functionType) with a spike-
            trian produced by the neuron (indicated here again in self.functionType). It is important to note that many
            times the indices in self.filterPositions are duplicated, this occurs when the functions (as indicated in
            self.functionType) are the identical except for the "timeDelay" value. For reasons of memory use, we try not
            to duplicate filterd spike-trains where not necessary. This is done by storing only one filtered spike train
            but then actually grabbing data from different points in time. The delay in time-steps is indicated in
            self.filteredStepsDelayed (see below). NOTE, if this value is -1, this means we do not reference the
            .filteredArrays and rather this is a bias term.
        self.filterStepsDelayed: the number of time-steps backwards (in condInenFuncFilteredManager.filteredArrays) we
            look backwards because of some function delay. Note that many functions that constitute the GLM can be
            viewed as a time-delayed Gaussian (for example). Instead of storing, what literally is the same filtered
            values over and over again in condInenFuncFilterManger.filterdArrays, we only store them once with a integer
            that indicates how many time-steps backwards we need to look. NOTE this value is the same as that stored in
            condInenFuncFilterManger.numbStepsToStore
        self.filterCoefficients: 2D numpy.array, Nx1, where N is the number of functions, this is how we store all the
            coefficients that need to get updated during training.
        self.dataMatrix: nXm 2D numpy.array, where the number of rows is equal to the number of data samples, and the
            number of columns equal to the number of coefficients. When this is created, it is ASSUMED,
            that the data samples come from consecutive time-steps, with the last data-point coming from
            the point in time last recorded here in self.spikeBoolTimeArray.
        self.spikeArray: for each time-step (as noted in self.spikeTimeArray), we have either a 0 or 1 to indicate
            no-spike or spike respectively
        self.spikeTimeArray: the associated time of the spike given in self.spikeBoolArray

        :param neuronGroup: integer, the neuron-group the neuron resides within
        :param neuronId: integer, the neuron index within the group
        :param dt: double, units: seconds, the time-step used for all the CIFs
        """
        self.neuronGroup        = neuronGroup
        self.neuronId           = neuronId
        self.dt                 = dt
        self.functionType       = [] # List of function types, for example [ ["gausLinSelf', a, b, c, timeLength, dt, timeDelay], ...]
        self.filterPositions    = [] # The list of indicies, where each index value points to a filterd neuron spike train stored in condInenFuncFilterManager.filteredArrays
        self.filterStepsDelayed = [] # The list of time-step delays used for each of the filtered data.
        self.filterCoefficients = numpy.array([[]]) # The list of coefficients that go along with each filter.
        # This following array will change every time we update data to be used to predict the states
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
        Add a new functionType/filter to the CIF object. NOTE, this method includes checking, thus before actually
        adding the new information, we make sure it is not all ready included. If it already is, then the given info is
        NOT added as no duplicates are allowed.

        When we say we are "adding" a filter/function, we are not here adding the actual array or numerical quantity or
        expression. Rather we are adding a textual description of the filter ('functionType'), and the index of
        the target neuron's fitlered spikes as stored in the condInenFuncFilterManager.filteredArrays ('filterPosition')
        , and lastly we are storing the 'stepsDelayed' which are the time-steps of delay that cannot be stored
        implicitly in the filter.

        :param functionType: list of lists for the following forms, where each indicates the functional form of one of
            the components in the GLM, the acceptable forms allowed are as follows:
            ['bias']
            ['gausLinSelf', a, b, c, timeLength, dt, timeDelay]
            ['gausLinSyn', otherNeuronGroup, otherNeuronId, a, b, c, timeLength, dt, timeDelay]
            ['gausQuadSynSelf', otherNeuronGroup, otherNeuronId, a, b, d, timeLength, dt, timeDelaySelf, timeDelayOther]
        :param filterPosition: integer, index to a row of filtered spike data that is stored in some
            condInenFucFilterManager object (see condInenFuncFilterManager.filteredArrays), wherein that data is
            needed to create the CIF. NOTE, if this value is -1, this means we do not reference the .filteredArrays and
            rather this is a bias term.
        :param stepsDelayed: integer, the number of time-steps backwards (in reference to filtered data stored in some
            condInenFuncFilterManager.filterdArrays) that we grab the data from for the current time-step.

        :return: N/A
        """
        if not self.doesFunctionExist(functionType):
            self.functionType.append(functionType)
            self.filterPositions.append(filterPosition)
            self.filterStepsDelayed.append(stepsDelayed)
            if self.filterCoefficients.size > 0:
                self.filterCoefficients = numpy.vstack((self.filterCoefficients, 0.))
            else:
                self.filterCoefficients = numpy.array([[0.]])


    def appendSpikeArray(self, spikeTimes, startTime, endTime, dt):
        """
        DESCRIPTION
        For each neuron that has a CIF, we keep track of its spikes in memory, new spikes are appended here. This may
        be a bit of a memory suck as this likly replicates much of the memory in condInenFuncFilterManager object,
        but it is here for ease of use and manipulation.

        Here both self.spikeArray and self.spikeTimeArray are appended. Data is only deleted if it overlaps currently
        stored values. We do not store spike-times explicity in this class, rather we have a boolean array (stored in
        self.spikeArray) and an array of equal size that indicates the time step (stored in self.spikeTimeArray)

        IMPORTANTLY, as in with the convolution functions, we do not store or include the 'endTime' time-step. That is
        all the spikes that occur from t=endTime-dt -to- t=endTime are stored in the time-step = endTime-dt (this is
        specifically the case when the duration of time given is an integer multiple of dt).

        NOTE: This method should be called when new simulation is added and filtered so that the data here is up to
        date.

        NOTE: methods provided in this class to handle cleanup and deletion of memory in order to reduce memory
        overhead.

        :param spikeTimes: 1D numpy.array, of neuron spike times
        :param startTime: double, units - seconds
        :param endTime: double, units - seconds

        :return: N/A
        """
        # Create boolean array in place of the spike times, remember we do not include the endTime time-step! NOTE:
        # rounding problems with numpy.arange have been problematic, so going to use numpy.linspace instead. The
        # arange function can lead to more pronouced rounding errors. Also we need to use numpy.round, before using
        # the int() command, as int() does not round as expected. When we make the time array we do not want to include
        # the 'endTime'.
        #timeArray = numpy.arange(start=startTime, stop=endTime, step=dt)
        timeArray = numpy.linspace(startTime, endTime-dt, int(numpy.round((endTime-startTime)/dt)))
        boolArray = numpy.zeros(timeArray.shape[0])
        indsArray = numpy.round((spikeTimes-startTime)/dt) # Use round() function before passing into int() function
        # Make sure if there is a spike that has a spike-time of end-time, we make it slightly less so it rounds down
        indsArray[indsArray >= boolArray.size] -= 1
        boolArray[indsArray.astype(int)] = 1.
        # Append/update self.spikeArray and self.spikeTimeArray
        if self.spikeArray.size > 0:
            # Find the next expected time-step that is thought to be seen occuring next in the input data. Depending on
            # on what this is will change how the new data is added to the old stored data.
            expectedNextTime = self.spikeTimeArray[-1] + dt
            # Check if there is missing time between the stored and given data, thus we have to add filler time and
            # filler data. Then append the new information. NOTE, because of rounding issues we cannot directly compare
            # floats in the if-statements conditional.
            if (timeArray[0] - expectedNextTime) > dt*.95:
                #fillerTimes = numpy.arange(start=expectedNextTime, stop=startTime, step=dt)
                fillerTimes = numpy.linspace(start=expectedNextTime, stop=startTime-dt, num=int(numpy.round((expectedNextTime-startTime)/dt)))
                fillerZeros = numpy.zeros(fillerTimes.size)
                # Because of rounding issues, we aagain want to double check and make sure we are updating as expected.
                if fillerZeros.size>0:
                    self.spikeArray = numpy.hstack((self.spikeArray, fillerZeros, boolArray))
                    self.spikeTimeArray = numpy.hstack((self.spikeTimeArray, fillerTimes, timeArray))
                else:
                    self.spikeArray = numpy.hstack((self.spikeArray, boolArray))
                    self.spikeTimeArray = numpy.hstack((self.spikeTimeArray, timeArray))
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


    def createFeasibleFilteredDataSet(self, filteredArrays):
        """
        DESCRIPTION
        Here the filtered-data is prepared to be used as data for prediction of the this neuron's spikes. As this class
        does not store that data (filtered spike trains from this and other neurons) it must be given here. There are
        likely more efficient ways to do this (not passing data into this method) but it suffices for now.

        We are parsing  'filteredData' which is (or should be) litteraly the list of numpy.array objects stored in
        condInenFuncFilterManager.filteredArrays.

        When we are creating this data, we are ASSUMING and that all the data in 'filteredArrays' is all up to date, and
        that all the most recent spike data has been used

        :param filteredArrays: list of 1D numpy.arrays, each numpy.arrray is some neuron's spikes filtered by some
            filter.

        :return: N/A
        """
        # We want to know how many time-steps that CAN be possibly be used for each of the filtered data arrays. In
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
                maxSteps[ind] = 10000000000
        maxStepsAllowed = int(numpy.amin(maxSteps)) # At most we take the minimum value.
        # Now given the maximum number of time-steps where we have a complete data set, we need to actually collect the
        # data and put into the correct form, in some N-D numpy.array([[]]) where each row is to be a sample of data,
        # and then fill the data matrix, this matrix of data will be self.dataMatrix
        if maxStepsAllowed > 0:
            self.dataMatrix = numpy.ones((maxStepsAllowed, len(self.filterPositions)))
            # Iterate over each column, which are the different filtered spike trains.
            for i in range(self.dataMatrix.shape[1]):
                # If we have a bias term, we just keep the ones, otherwise we put the correct data in the arrays. Data
                # will be grabed from the "filteredArrays" input (list of numpy.array).
                if self.filterStepsDelayed[i]>0:
                    self.dataMatrix[:, i] = filteredArrays[self.filterPositions[i]][-(maxStepsAllowed+self.filterStepsDelayed[i]):-self.filterStepsDelayed[i]]
                elif self.filterStepsDelayed[i]==0:
                    self.dataMatrix[:, i] = filteredArrays[self.filterPositions[i]][-maxStepsAllowed:]
            # Now flip the matrix so that we have the lastest time along the bottom row, and the data from the earlier
            # times along the first row.
            # self.dataMatrix = numpy.flipud(self.dataMatrix)
                # AJL July 2016, this has been changed because the latest, most recent times are already along the
                # bottom row
        else:
            # If there is no data, we just leave an empty matrix.
            self.dataMatrix = numpy.array([])


    def generateSpikeProbabilities(self):
        """
        DESCRIPTION
        Here the probabilities of this neuron's spike(s) is calculated. This calculation is based on data stored in
        self.dataMatrix which must have already be calculated. The data in self.dataMatrix is organized such that each
        row is data input for one time-step and each column represents a different component of the CIF.

        It is assumed that self.dataMatrix has been updated and given the most up to date data.

        :return: the probability of this neuron spiking at each time-step (number of time-steps =
        self.dataMatrix.shape[0] which is the given number of data points available)
        """
        # Numpy is weird, they use the method 'dot' to handle matrix multiplication.... ugg I do not like that
        if self.dataMatrix.shape[0]>0 and self.dataMatrix.shape[1]>0:
            glmOutput = numpy.dot(self.dataMatrix, self.filterCoefficients)
            glmOutput = 1./(1. + numpy.exp(-1.*glmOutput))
        else:
            glmOutput = numpy.array([])
        return glmOutput


    def updateCifCoeffsSGD(self, learningRate):
        """
        DESCRIPTION
        Assuming self.dataMatrix, self.spikeArray, and self.spikeTimeArray are all up to date (self.dataMatrix gets
        updated in a different method that the other two) we go ahead and run stochastic gradient descent to update the
        CIF coefficients in the GLM.

        IMPORTANT, it is ASSUMED that the data stored in self.dataMatrix, self.spikeArray, and self.spikeTimeArray all
        have been updated with the same data. This should be the case if the input data from the simulation is used to
        (1) update the filter-manager class obj (condInenFuncFilterManager.filterSpikeTrainData(), (2) then using the
        same data update self.spikeArray and self.spikeTimeArray using self.appendSpikeArray(), and then subsequently
        calling self.createFeasibleFilteredDataSet(). These actions should be automated in the
        condItenFuncManger.parseNewSpikeTrainData() method.

        :param learningRate: double, controls rate at which SGD updates

        :return: N/A
        """
        t0 = time.clock()
        self.filterCoefficients = updateCifCoeffsSGD.updateCifCoeffsSGDcy(self.dataMatrix, self.spikeArray,
                                                                             self.filterCoefficients, learningRate)
        tfirst = time.clock() - t0
        myStr = "Time to CIF update = " + str(tfirst)
        print(myStr)


    def clearSpikeArray(self):
        """
        DESCRIPTION
        This is meant to clear out spike memory of the neuron in question. This should be typically run after
        probabilities are calculated and after CIF coefficients are updated. Then we will not need the stored data
        any longer

        :return: N/A
        """
        self.spikeArray = numpy.array([])
        self.spikeTimeArray = numpy.array([])