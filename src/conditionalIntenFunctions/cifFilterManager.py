# A. Lonsberry
# September 2016 - Taken from sandbox (condIntenFunc_0_3.py)
#
# DESCRIPTION
#   The class(es) below are designed to handle storing filters (that comprise the individual CIFs), parsing spike trains
#   filtering the spikes, and storing the results.


import numpy
import matplotlib.pyplot as plt
import numpy.matlib
import cifUtilityFunctions

class condInenFuncFilterManager:
    """
    DESCRIPTION
    Class: manages different filter-functions, filters spike-trains, and stores results.

    General Guidelines on usage:

    (1) Adding new filters: in order to add a new filter/function for some neuron, one will typically first call
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
        self.filterIdentities: list of list, [ [filter 1 info], [filter 2 info], ... ], where in each sub-list the
            following variants are acceptable:
            ['gaussian', a, b, c, timeLength, dt, timeDelay, truncation, filterIndex] where 'filterIndex' refers to the
            index or position of the actual filter as located within the self.filters list.
        self.filters: list of 1D numpy.array, filters where each filter is something we convolve the spike-trains that
            have been recorded.
        self.filtersDt: list of time-step values for each filter, this list should be as long as the self.filters list.
        self.filtersBackUpTime: list of the amount of time we want to include backwards in order to make the filters
            accurate. This inclusion of spikes over some backwards horizon is necessary in order insure the accuracy
            of the output for the interval of time we actually care about.
        self.neuronToFilter: 2D numpy.array, an example row being [ [neuron-group, neuron Id, filterIndex], ... ]
        self.numbStepsToStore: 2D numpy.array, the number of time-steps we need to keep for the next time we need to
            run calculations, essentially this is the maximum delay we need (July 2016, used in garbage collection and
            debugging)
        self.actualDelay: 2D numpy.array, tells us how much delay each filter has (July 2016, I think this is now only
            for debugging)
        self.filteredArrays: a list of 1D numpy.array, where in each array has the filtered spike trian as indicated in
            self.neuronToFilter
        self.filteredArraysTimeValues: list of 1D numpy.array, where in each array has the time-step of the filterd
            value in self.filteredArrays. Thus this list and self.filteredArrays are essential to each other. One gives
            the actual fitlered value, and the other the time in which that value occurs at.
        self.spikeTrainMemory:  a list of 1D numpy.array, where in each array we have all the spike times, [ [t1, t2, ...
                               tn], ... ]. The old spike times that are kept around.
        """
        # The size of the following list is based on the number of unique fiters and nuron combos
        self.filterIdentities     = [] # Association between filter-identity/filter-type to an actual filter in self.filters
        # The size of the following lists is eqaul to the number of unique convolution filters needed
        self.filters              = [] # The actual convoltuion filters which are 1D numpy.arrays
        self.filtersDt            = [] # The time-step associate to filter in self.filters
        self.filtersBackUpTime    = [] # The amount of time before we need to include backwards for each filter in self.filters
        # The size of the following are equal to the number of unique neuron-convoFilter pairs
        self.neuronToFilter            = numpy.array([[]]) # The association between filters and neurons, numpy.array([ [neuron-group, neuron Id, filterIndex], ...]
        self.numbStepsToStore          = numpy.array([[]]) # The number of time-steps we need to keep for the next time we run calculations, essentially this is the maximum delay we need
        self.actualDelay               = numpy.array([[]]) # Realy just kept as a measure to debug, this tells us how much of a delay each filter has.
        self.filteredArrays            = [] # The filtered spike-train output at each time-step resultant from convolution operation
        self.filteredArraysTimeValues  = [] # The time-step corresponding to each filtered output time-step, should be the same size as the filterArray
        self.spikeTrainMemory          = [] # The old spike times we want to keep for each neuron


    def getFilterIndex(self, filterInfo, timeDelay=None, truncation=None):
        """
        DESCRIPTION
        Given a filter info, along with either a 'timeDelay' or a 'truncation' we will find if the filter already exists
        and should it exist, this will return the position of the associated filter (a 1D numpy.array) within the list
        of said filters denoted as self.filters.
        :param filterInfo: a list of filter information, where the following variants are acceptable,
            ['gaussian', a, b, c, timeLength, dt]

        :return: filter index (integer) -or- -1
        """
        toReturn = -1 # Initialize the output
        if filterInfo[0]=='gaussian':
            for aList in self.filterIdentities:
                if aList[0:6]==filterInfo:
                    if timeDelay!=None:
                        if aList[6]==timeDelay:
                            toReturn = aList[-1] # We always store the filterIndex last
                            break
                    else:
                        if aList[7]==truncation:
                            toReturn = aList[-1] # We always store the filterIndex last
                            break
        return toReturn


    def getFilterTruncation(self, filterInfo, timeDelay=None, filterIndex=None):
        """
        DESCRIPTION
        Given a filter info, along with either a 'timeDelay' or a 'filterIndex' we will find if the filter already exists
        and should it exist, this will return the truncation associated with the filter.

        :param filterInfo: a list of filter information, where the following varients are acceptable
            ['gaussian', a, b, c, timeLength, dt]

        :return: filter truncation -or- -1
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
        function-type describes the filter stored in the list of such filters, self.filters[fitlerIndex].

        :param filterIndex: integer, the index of a particular filter in self.filters, (self.filters is a list of 1D
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


    def getNeuronAndFilterIndex(self, neuronGroupId, neuronId, filterIndex):
        """
        DESCRIPTION
        Given the the quantities to identify the neuron, which are 'neuronGroupId' and 'neuronId', and along with the
        'filterIndex' which corresponds to the location of some filter within the list of filters, self.filters, we
        find if this combination already exists and if it does what is the index within self.neuronToFilter.

        :param neuronGroupId: integer, the neuron-group ID
        :param neuronId: integer, the neuron ID within its group
        :param filterIndex: integer, the index of the filter as found in self.filters

        :return: index of the combination in self.neuronToFilter -or- -1 should it not exist
        """
        inds = numpy.where((self.neuronToFilter == (neuronGroupId, neuronId, filterIndex)).all(axis=1))
        if inds[0].size > 0:
            dummyOut = inds[0][0]
        else:
            dummyOut = -1
        return dummyOut


    def addFilter(self, filterInfo):
        """
        DESCRIPTION
        Adds a new filter, where a "filter" here is something specifically that is will be convolved with spike trains.
        The full filtered is described by the input 'filterInto'. Importantly, this method includes internal
        checks to make sure no duplicates are added. Thus, if a filter is specified in the input that is already stored,
        a duplicate will NOT be added.

        This function DOES NOT handle matching of neurons to filters. It is simply here to add new filters when needed
        and which not already represented. To see how neurons are matched with filters see the methods below,
        self.addNeuronAndFilter() and self.initFilterRecorders().

        An example of 'filterInfo is ['gaussian', a, b, c, timeLength, dt, timeDelay]. Now remember, there maybe more
        self.filterIdentities than actual self.fitlers. This is because many of those filters described in the
        self.filterIdentities are actually produce the same output just delayed by different time-steps (The filters
        have the same shape)!

        Some quick notes on time-delays, time-lengths, and truncations. A filter (as of the time of this writing July
        2016) has a timeLength which is 1/2 of the filter length. If we assume the filter is centered at some origin,
        the 'timeLength' is from the origin to the most positive value in the domain. The time-delay is the amount of
        time the filter is delayed. That is how far the origin of the filter is shifted to the left. Truncation is not
        given, it is calculated below as min(timeLength, timeDelay) and indicates the actual domain of the filter which
        is [-timeLength, truncation].

        IMPORTANTLY if the truncation value is > 0, then a time-delay is implicitly store in the filter. If the
        time-delay is greater than the time-length, then we will also internally create 'timeStepsDelayed' which is to
        store the amount of delay greater than that stored implicitly (time-length).

        :param filterInfo: a list of filter information, where the following variants are acceptable,
            ['gaussian', a, b, c, timeLength, dt, timeDelay]

        :return: list with the values, [truncation, filterIndex, timeStepsDelayed]
        """
        if filterInfo[0]=='gaussian':
            # First determine the number of time-steps we need to look backwards in memory of the filtered spike-trains.
            # That is some filters are time-delayed Gaussians. Each of the filters created implicitly have a delay
            # emmbedded. However if the the the time-delay is longer than the length of the filter ("timeLength") then
            # then there will be an addition of time-steps to look backwards.
            timeDelay  = filterInfo[6]
            timeLength = filterInfo[4]
            truncation = min(timeLength, timeDelay)
            timeStepsDelayed = int((timeDelay - truncation)/filterInfo[5]) # This is the number of time-steps backwards we need to look backwrds.
            if self.getFilterIndex(filterInfo[0:6], timeDelay=filterInfo[6]) != -1:
                # This filter already exists, which means we only have to grab the relevant output details which in this
                # case is the size of the filter's truncation, and the filter-index.
                dummyList = [self.getFilterTruncation(filterInfo[0:6], timeDelay=filterInfo[6]), self.getFilterIndex(filterInfo[0:6], timeDelay=filterInfo[6]), timeStepsDelayed]
            else:
                # Filter of that type does not exist, Remember, though some filters are different,
                # they essentially are the same filter convolved with a spike train, but with differnet delay values.
                if self.getFilterIndex(filterInfo[0:6], truncation=truncation) != -1:
                    # Filter exists with different time-delay (but really same filter/convolution operation). Record the
                    # filterInfo, but do not create a new filter for convolution.
                    newFilterIndex = self.getFilterIndex(filterInfo[0:6], truncation=truncation)
                    filterInfo.append(truncation) # Append to the filter information the truncation time
                    filterInfo.append(newFilterIndex) # Append to the filter the index of the filter in self.filters
                    self.filterIdentities.append(filterInfo)
                else:
                    # In this case we have to create a new filter because it does not already exist. Firstly we want
                    # to determine the filter-backup-time. This is the amount of time over some backwards horizon
                    # wherein we need to store neuron spikes in memory so the filter values over the time we care
                    # about is correct!
                    filBackUpTime = truncation - (-1.*timeLength)
                    self.filters.append(cifUtilityFunctions.truncatedGaussianFilterGenerator(filterInfo[1], filterInfo[2], filterInfo[3], filterInfo[4], filterInfo[5], truncation))
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
        This class contains/stores filters (which are things we convolve spike trains with). It also stores the results
        of convoltuion between said filter with neuron spike-trains. The question that remains is what filters need to
        be applied to which spike trains? This problem is solved with self.neuronsToFilter which is a 2D numpy.array of
        size NX3, where each row has the format [neuron-group, neuron id, filterIndex]. The filter index in this case
        is the index of the filter stored here in self.filters.

        Here a new row is added to self.neuronsToFilter only if, that row does NOT already exist. That is this method
        will internally check to make sure this is not a duplicate entry. If it is, then a new row will not be added.

        Other than self.neuronsToFilter, self.numbStepsToStore and self.actualyDelay will be updated as well though
        those two variables (as of July 2016) seem to be more for debugging purposes.

        :param neuronGroupId: integer, the neurons's group id
        :param neuronId: interger, the neurons index in its group
        :param filterIndex: the index in self.neuron
        :param timeStepsDelayed: how many steps backwards in filtered data we have to look, some filters have implicit
            delays and will not actually have any steps to look backwards.

        :return: integer, index of self.filteredArrays where the filterd spike trains will be stored.
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
        Once all the filters needed have been added via self.addFilter(), and all the neurons that need
        their spike trains filtered have been indicates via the class method self.addNeuronAndFilter(), then this
        function MUST BE CALLED.

        This function initializes the class variables that store all the spikes, and filtered spike trains. In
        particular this will be self.filteredArrays, self.filteredArraysTimeValues, and self.spikeTrainMemory. The
        first of these two are lists of 1D numpy.array objects that store the filtered and time of filtered values
        respectively. The last class variable will store again a list of numpy.array objects, but each numpy.array will
        sotre spike times of the appropriate neuron as indicated by self.neuronToFilter. Thus self.neuronToFilter
        corresponds to each of this initialized variables.

        It is only allowed to call this operator once. Every time it gets called if will overwrite all formally stored
        data.

        NOTE: some of the memory will be duplicates (at the moment) where we are repeatedly storing the spikes from the
        same neuron over and over again.

        :return: N/A
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

        IMPORTANT to note, we assume data is given each time for ALL neurons that are being tracked in
        self.neuronsToFilter! This is very important. If there is no spike data for one of those neurons, it is assumed
        that said neuron produced zero spikes from startTime to endTime.

        IMPORTANT to note, the last time-step at time=endTime is NOT included or stored. This is because of the way we
        time-step through data. Each time-step includes all the spikes up to, but not including the starting point of
        the next time-step.

        IMPORTANT to note, the brain2.SpikeMonitor objects must index their neurons, in the same way that the neurons
        have been index in all the filter-manager variables. For example ,if a neuron belongs to neuron-group 2, with
        neuron-id of 5, then we will expect the brain2.SpikeMonitor object for neuron-group2 to have the spikes stored
        in said object ot be associated with index 5.

        :param listOfSpikeObjs: list of brian2.SpikeMonitor class objects
        :param listOfNeuronGroups: list of neuron-groups, simply associating each SpikeMonitor object with a particular
            neuron-group and the neurons in said group
        :param startTime: double, the start-time of the added data
        :param endTime: double, the end-time of the added data

        :return: N/A
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
            if nGrp in listOfNeuronGroups: # Check if the input data has the neuron group represented or not
                spkObjInd = listOfNeuronGroups.index(nGrp) # Indicates which of the brian2.SpikeMonitor objects to grab from 'listOfSpikeObjs' list
                spkObj = listOfSpikeObjs[spkObjInd] # Grab the correct brian2.SpikeMonitor object
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
                filteredSpikeTrainList = cifUtilityFunctions.convolveSpikeTimesWthBackUp(self.filters[filt], spikeTimes, startTime, endTime, self.filtersDt[filt], backupTime=self.filtersBackUpTime[filt])
                # Call functions to append filtered data, and spike data to memory.
                self.appendFilteredArrays(ind, filteredSpikeTrainList[0], filteredSpikeTrainList[1])
                self.appendSpikeTrainMemory(ind, spikeTimes)


    def appendFilteredArrays(self, ind, filteredArray, filteredArrayTimes):
        """
        DESCRIPTION
        A method, called generally when given new spike data that has been filtered, to append self.filteredArrays
        and self.filteredArraysTimeValues.

        IMPORTANT to note, data is deleted if there is an overlap in time-arrays.

        IMPORTANT to note, if there is a gap in time between the stored memory and the given data, that time is filled
        in with blank data to make there is no gap in time.

        :param ind: integer, the index of the current neuron / filter combination stored in self.neuronToFilter where
        :param filteredArray: numpy.array([]), the filtered data
        :param filteredArrayTimes: numpy.array([]), associated time of filteredArray
        :return: N/A
        """
        if filteredArray.size>0: # Check if the input data is valid
            #Find the time-step!
            dt = self.filtersDt[self.neuronToFilter[ind, 2]]
            # Check if the stored array previously has values.
            if self.filteredArraysTimeValues[ind].size>0:
                # Find the next expected time-step that is thought to be seen occuring next in the input data. Depending
                # on what this is will change how the appending process works.
                expectedNextTime = self.filteredArraysTimeValues[ind][-1] + dt
                # Check if there is missing time between the stored and given data, thus we have to add filler time and
                # filler filtered data. Then append the correct self.filtereArrays and self.filterArraysTimeValues. Note
                # the code below may seem a bit odd, but because of rounding errors, when numpy.arange() is used (it is
                # used to create 'filterArrayTimes') rounding errors tend to be involved so one cannot expect these
                # to be exactly correct. Thus we introduce numpy.linspace instead.
                # if filteredArrayTimes[0]>expectedNextTime:
                if filteredArrayTimes[0]-expectedNextTime > dt*.95:
                    #fillerTimes = numpy.arange(start=expectedNextTime, stop=self.filteredArraysTimeValues[ind][0], step=dt)
                    fillerTimes = numpy.linspace(start=expectedNextTime, stop=filteredArrayTimes[0]-dt, num=int(numpy.round((filteredArrayTimes[0]-expectedNextTime)/dt)))
                    fillerZeros = numpy.zeros((1, fillerTimes.size))
                    # Again we need to double check for rounding issues once more.
                    if fillerZeros.size>0:
                        self.filteredArrays[ind] = numpy.hstack((self.filteredArrays[ind], fillerZeros, filteredArray)) #Filtered spike train
                        self.filteredArraysTimeValues[ind] = numpy.hstack((self.filteredArraysTimeValues[ind], fillerTimes, filteredArrayTimes)) #Associated time arrays
                    else:
                        self.filteredArrays[ind] = numpy.hstack((self.filteredArrays[ind], filteredArray))  # Filtered spike train
                        self.filteredArraysTimeValues[ind] = numpy.hstack((self.filteredArraysTimeValues[ind], filteredArrayTimes))  # Associated time arrays
                # Check if there is an overlap in time, between stored data and that of the newly given data. This is
                # not desireable, best is that given data is a continuation of some simulation. However, if this is the
                # case, we clear out part of the memory in question, and replace it. It is UP TO THE USER to make sure
                # that this operation makes sense. Note the code below may seem a bit odd, but because of rounding
                # errors, when numpy.arange() is used (it is used to create 'filterArrayTimes') rounding errors tend to
                # be involved so one cannot expect these to be exactly correct.
                # elif filteredArrayTimes[0]<expectedNextTime:
                elif expectedNextTime-filteredArrayTimes[0] > dt*.05:
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
                # appending process. This is the most desirable situation where new data given is immediatly follows
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
        If this method is called, it is assumed that all the CIF prediction and updating for some epoch has concluded or
        is not needed. This method will go through and shrink the numpy.arrays in self.filteredArrays and in
        self.filteredArraysTimeValues to the smallest sizes possible as denoted by the values stored in
        self.numbStepsToStore.

        :return: N/A
        '''
        # Itearte over all the neuron / filter combinations saved in self.neuronToFilter
        for ind in range(self.numbStepsToStore.size):
            dStps = self.numbStepsToStore[ind] # The number of time-steps we want to save in the filtered arrays.
            if dStps > 0:
                self.filteredArrays[ind] = self.filteredArrays[ind][-dStps:]
                self.filteredArraysTimeValues[ind] = self.filteredArraysTimeValues[ind][-dStps:]
            else:
                self.filteredArrays[ind] = numpy.array([])
                self.filteredArraysTimeValues[ind] = numpy.array([])


    def clearFilteredArrays(self, ind, startTime=0., stopTime=-1):
        """
        DESCRIPTION
        A method to clear some self.filtered arrays, for a particular neuron / filter indicated by the postition in
        self.neuronToFilter by the input 'ind', over some time-range given by starTime - infinity, or startTime to
        stopTime.

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

        :param ind: the index of the neuron-group, neuron, neuron-filter combination as stored in self.neuronToFilter
        :param spikeTimes: numpy.array([]), the set of spike times to append
        :return: N/A
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
        If this method is called, it is assumed that all the CIF prediction and updating for some epoch has concluded or
        is no longer needed. This method will go through and shrink the numpy.arrays in self.spikeTrainMemory to the
        smallest sizes possible as denoted by the value stored in self.filtersBackUpTime, that is each filter
        (self.filters) has a filter-backup time (self.filterBackupTime) that denotes this amount of time.

        This function will ONLY keep the data from [endTime-filterBackupTime, endTime].

        :return: N/A
        '''
        for ind, filtInd in enumerate(self.neuronToFilter[:, 2]):
            self.spikeTrainMemory[ind] = self.spikeTrainMemory[ind][self.spikeTrainMemory[ind] >= endTime-self.filtersBackUpTime[filtInd]]
            self.spikeTrainMemory[ind] = self.spikeTrainMemory[ind][self.spikeTrainMemory[ind] <= endTime]


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

        :return: N/A
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


    def clearAllMemory(self):
        """
        DESCRIPTION
        Here we simply clear out all the memory and restart
        :return:
        """
        self.initFilterRecorders()