# A.Lonsberry
# October 2016
#
# DESCRIPTION
#   Here I am going to have methods to create and save data for this experiment. The methods below are VERY hardcoded,
#   but rely on other methods that are generic. They are hardcoded so we can generate very specific data (original data,
#   decomposed original data, and then spiking from analog)

import matplotlib.pyplot as plt
import numpy

# I had a difficult time getting things to import, so now I have however, go about fixing this by going up a directoy
# which I think I can do in the following way.
import os
import sys
import pickle
sys.path.append(os.path.abspath("../"))

# Call the python file that has the 1D mass simulator in it.
from dataGenerators import movingPointMass1D
from spikeDataGenerators import decompDimensions
from spikeDataGenerators import analogToSpikes

########################################################################################################################
# METHODS FOR CREATING 1D MOVING MASS DATA (SMOOTH, DECOMPOSED, SPIKING)
########################################################################################################################

def genAndSaveMoving1DMassData(saveName='movingPointMassData/pointMassData000.pkl',Iterations=10):
    """
    DESCRIPTION
    Here we call the methods repeatedly to make instances of the 1D moving mass. Of these n examples, we then put them
    into a list, and then save them in the appropriate folder.
    :return: N/A (we save data to a file), the list saved is [listOfTrials, xmin, xmax, vmin, vmax, amin, amax, dt, tmax]
    """
    #How many iterations we want to include
    #Iterations = 10 # No more hard coding! This is now an input into the function.
    xmin = 0.0
    xmax = 20.0
    vmin = 1.0
    vmax = 5.0
    amin = 1.0
    amax = 2.0
    dt   = .001
    tmax = 20.0
    dataOut = []
    for i in range(Iterations):
        funcOuts = movingPointMass1D.randomMassMovement(xmin, xmax, vmin, vmax, amin, amax, dt, tmax)
        dataOut.append(funcOuts)
    toSave = [dataOut, xmin, xmax, vmin, vmax, amin, amax, dt, tmax]
    outputFile = open(saveName, "wb")
    pickle.dump(toSave,outputFile)
    outputFile.close()

def loadAndPlot1DMassData(dataFile='movingPointMassData/pointMassData000.pkl'):
    """
    DESCRIPTION
    We load in some of the saved data from 1D moving mass, and then plot what those trajectories look like. This method
    is given primarily as a means of debugging and checking.
    :param dataFile: The file where the 1D mass movement data is stored.
    :return: N/A we just plot stuff out.
    """
    # Load the data back
    inputDataFile = open("movingPointMassData/pointMassData000.pkl", "rb")
    dataOut = pickle.load(inputDataFile)
    inputDataFile.close()
    # Iterate over the different saved trajectores and plot out the results.
    for i in range(len(dataOut[0])):
        plt.figure(i)
        plt.plot(dataOut[0][i][1],dataOut[0][i][0])
    plt.show()


def decomposeMoving1DMassData(dataFile='movingPointMassData/pointMassData000.pkl', saveName='movingPointMassData/pointMassDataDecmp000.pkl' ):
    """
    DESCRIPTION
    This is made to take 1D point mass trails, given the file name is given, and then subseqently make normal functions
    to decompose the original 1D doamin, and then actually decompose the domain, and finally save the results.
    :param dataFile: string, the location where we have the data stored
    :param saveName: string, the location where we want to save the results
    :return: N/A (we save data to a file), the list saved is [list of segmented 1D trails (numpy.arrays), gCenters (
            numpy.array), b (scalar), file name of original data (string)]
    """
    # Load in the saved 1D moving mass. Inside the loaded list of lists, we will have a list that has a bunch of
    # seperate 1D numpy.arrays that represent the position of the mass through different random trials.
    inputFile = open(dataFile, "rb")
    trackingData = pickle.load(inputFile) # data in list like [dataOut, xmin, xmax, vmin, vmax, amin, amax, dt, tmax] see genAndSaveMoving1DMassData()
    inputFile.close()
    trackedArrays = trackingData[0] # This is a list of lists where in each we have [ [positionArray, timeArray], [positionArray, timeArray], .... [positionArray, timeArray] ]
    print(len(trackedArrays))

    # Now we need to generate the Gaussian functions that will be used to decompose the 1D data into separate parts
    # that will represent inputs from individual neurons.
    xStart = 0.0
    xStop = 5.0
    nGaus = 10
    bInitial = .05
    dt = 0.001
    cfds = decompDimensions.minVarOverlappingGaussian(xStart, xStop, nGaus, bInitial, dt)
    gCenters = cfds[0]
    b = cfds[1]

    # Now we need to parse the loaded data with the generated Gaussians
    #seperateInputValuesAlongGaussians(trackingData, gcenters, b, normalize=False)
    segmentedTrials = [] #This is where we will put the decomposed 1D mass trials.
    for i in range(len(trackedArrays)):
        decomped = decompDimensions.seperateInputValuesAlongGaussians(trackedArrays[i][0], gCenters, b)
        segmentedTrials.append(decomped)

    # Now we want to save the decomposed 1D arrays and the Gaussians that made them, and also cite the original
    # data that was used.
    toSave = [segmentedTrials, gCenters, b, dataFile, xStart, xStop, nGaus, bInitial, dt]
    outputFile = open(saveName, "wb")
    pickle.dump(toSave, outputFile)
    outputFile.close()


def loadAndPlotDecomp1DMassData(dataFile='movingPointMassData/pointMassDataDecmp000.pkl'):
    """
    DESCRIPTION
    We load in some of the saved data from 1D moving mass that have been decomposed into different receptive regions.
    Where we plot these to make sure they make sense. This is really only useful for debugging. NOTE: Make sure that the
    two files you are loading were generated by the same original data! (decomped applied to the correct randomly
    generated 1D mass movement data)
    :param dataFile: The file where the 1D mass movement data is stored.
    :return: N/A we are plotting here.
    """

    # Load in modules to handle the 3D plot (which I still do not well understand)
    from matplotlib.collections import PolyCollection as pc
    from mpl_toolkits.mplot3d import Axes3D

    # Load the data back (this is the decomposed version of the 1D moving mass data)
    inputDataFile = open(dataFile, "rb")
    dataOut = pickle.load(inputDataFile)
    inputDataFile.close()
    gCenters = dataOut[1] # The centers of the Gaussaians

    # Load in the original data (the filename is included in the loaded bit) with is the original 1D analog signal
    inputDataFile = open(dataOut[3], "rb")
    dataOrig = pickle.load(inputDataFile) # The original 1D mass movement data
    inputDataFile.close()

    # Now I need to plot these things out, iterate over the original 1D mass data.
    for i in range(len(dataOrig[0])):

        # Plot out the original data
        plt.figure(1)
        plt.plot(dataOrig[0][i][1], dataOrig[0][i][0])

        # Now plot out the decoped bits
        segmentedValues = dataOut[0][i]
        fig = plt.figure(2)
        ax = Axes3D(fig)  # Because I am using older version
        verts = []
        for j in range(dataOut[1].size):
            segmentedValues[0, j] = 0
            segmentedValues[-1, j] = 0
            # print(list(zip(segmentedValues[:,i],dArray)))
            verts.append(list(zip(segmentedValues[:, j], dataOrig[0][i][1])))
        poly = pc(verts)
        ax.add_collection3d(poly, gCenters, zdir='y')
        ax.set_xlim3d(0, 1.2)
        ax.set_zlim3d(0, 5)
        ax.set_ylim3d(0, 6)
        plt.show()


def decompedToSpikes1DMassData(dataFile='movingPointMassData/pointMassDataDecmp000.pkl', saveName='movingPointMassData/pointMassDataDecmpSpikes000.pkl'):
    """
    DESCRIPTION
    We load in the saved data, the decomposition of 1D moving mass data over n receptive feilds, and then convert each
    receptive field output into a set of spike trains.
    :param dataFile: string, the location where we have the decomp. data from the 1D movements stored.
    :param saveName: string, the location where we want to store the spiking version of the decomped data.
    :return: N/A we save the results
    """

    # I need to set the mechanism of spike generation! We could do it in a number of ways. For example by time-rescaling
    # For time-rescaling I am relying on the algorithms presented by Brown et al. "The TimeRescaling
    # Theorem and Its Application to Neural Spike Train Analysis" (2001) and their simple algorithm to make spikes. NOTE,
    # this is meant to handle making data look like it is generated from Poisson process, so there is a measure of
    # stochasticity involved here. Use another function if this is not wanted. Alternatively I have just linearly make
    # spikes based on inegrating each analog signal.
    #
    # I need to fill in the data here below, (Hardcoding) in order to determing what type of thing I am doing.
    #
    # Fill in spikeGenType: tells us which spike generation optiont to call. The options are "linearIntegral" or
    # "timeRescale"
    spikeGenType = "linearIntegral"
    analgSingnalScaling = 20.
    # Fill in the

    # Load the data back (this is the decomposed version of the 1D moving mass data)
    inputDataFile = open(dataFile, "rb")
    dataOut = pickle.load(inputDataFile)
    inputDataFile.close()
    segmentedTrialsList = dataOut[0] # list of numpy.arrays(time-steps X numb. decomps) of the decomposed 1D signal
    gCenters = dataOut[1]  # The centers of the Gaussaians
    b = dataOut[2]
    origFileName = dataOut[3]

    # I need to load the original filename to get the value of 'dt'
    inputDataFile = open(origFileName, "rb")
    origDataOut = pickle.load(inputDataFile)
    inputDataFile.close()
    dt = origDataOut[7]

    # Now make spikes
    segmentedSpikesList = []
    if spikeGenType=="linearIntegral":
        for i in range(len(segmentedTrialsList)):
            segmentedSpikesList.append(analogToSpikes.genSpikesLinearly(segmentedTrialsList[i], dt=dt, scaling=analgSingnalScaling))

    # spikeTrains = analogToSpikes.genSpikesWithTimeRescaling(segmentedValues, dt=.001, scaling=10.)
    # spikeTrains = analogToSpikes.genSpikesLinearly(segmentedValues, dt=.001, scaling=20.)

    # Save
    outputList = [segmentedSpikesList, dataFile, spikeGenType, analgSingnalScaling] # We also want to store the location of the originating decomped data
    outputFile = open(saveName, "wb")
    pickle.dump(outputList, outputFile)
    outputFile.close()


def loadAndPlotDecomp1DMassData(dataFile='movingPointMassData/pointMassDataDecmpSpikes000.pkl'):


    # Load in modules to handle the 3D plot (which I still do not well understand)
    from matplotlib.collections import PolyCollection as pc
    from mpl_toolkits.mplot3d import Axes3D

    #Load the data back (this is the spike version of the decomped values)
    inputDataFile = open(dataFile, "rb")
    dataOut = pickle.load(inputDataFile) # The saved list from turning decmped signals into spikes, [segmentedSpikeList, dataFileName]
    inputDataFile.close()
    segmentedSpikesList = dataOut[0]
    segmentedValuesFileName = dataOut[1]

    # Load the data back (this is the decomposed version of the 1D moving mass data)
    inputDataFile = open(segmentedValuesFileName, "rb")
    dataOut = pickle.load(inputDataFile) # The saved list from decmposing 1D mass movements [segmentedTrialsList, gCenters, b, dataFileName]
    inputDataFile.close()
    segmentedValuesList = dataOut[0]
    gCenters = dataOut[1]  # The centers of the Gaussaians
    origDataFileName = dataOut[3] # The save list of 1D mass movements

    # Load in the original data (the filename is included in the loaded bit) with is the original 1D analog signal
    inputDataFile = open(origDataFileName, "rb")
    dataOut = pickle.load(inputDataFile)  # The original 1D mass movement data [TrailsList, xmin, xmax, vmin, vmax, amin, amax, dt, tmax]
    inputDataFile.close()
    massMovementsList = dataOut[0] # The list of original mass movements

    # Now I need to plot these things out, iterate over the original 1D mass data.
    for i in range(len(massMovementsList)):

        # Plot out the original data
        plt.figure(1)
        plt.plot(massMovementsList[i][1], massMovementsList[i][0])

        # Now plot out the decomposed mass movements
        segmentedValues = segmentedValuesList[i]
        print(segmentedValues.shape)
        fig = plt.figure(2)
        ax = Axes3D(fig)  # Because I am using older version
        verts = []
        for j in range(segmentedValues.shape[1]):
            segmentedValues[0, j] = 0
            segmentedValues[-1, j] = 0
            # print(list(zip(segmentedValues[:,i],dArray)))
            verts.append(list(zip(segmentedValues[:, j], massMovementsList[i][1])))
        poly = pc(verts)
        ax.add_collection3d(poly, gCenters, zdir='y')
        ax.set_xlim3d(0, 1.2)
        ax.set_zlim3d(0, 5)
        ax.set_ylim3d(0, 6)
        #plt.show()

        #Now add spikes that need to be plotted on top as well which will be indicated with black lines
        verts = []
        segmentedSpikes = segmentedSpikesList[i]
        for j in range(len(segmentedSpikes)):
            if len(segmentedSpikes[j])>0:
                # Convert the list to a numpy array so we can plot stuff
                aSpikeTrain = numpy.asarray(segmentedSpikes[j])
                #
                x = numpy.zeros(aSpikeTrain.size*3)
                y = numpy.zeros(aSpikeTrain.size*3)
                ind = 0
                for k in range(aSpikeTrain.size):
                    x[ind] = aSpikeTrain[k]
                    y[ind] = 0.
                    ind += 1
                    x[ind] = aSpikeTrain[k]
                    y[ind] = 1.
                    ind += 1
                    x[ind] = aSpikeTrain[k]
                    y[ind] = 0.
                    ind += 1
                verts.append(list(zip(y,x)))
            else:
                # Need to give something to those without spikes or graph does not work out
                x = numpy.array([0., .001])
                y = numpy.array([.001, .001])
                verts.append(list(zip(y, x)))
        poly = pc(verts)
        ax.add_collection3d(poly, gCenters, zdir='y')
        plt.show()


def createDataDescriptionTxtFile(pMassFile=[], pMassDcmpFile=[], pMassDcmpSpkesFile=[]):
    """
    DESCRIPTION
    Here I want to create a .txt file that will describe the data saved in a particular set of files.
    :param pMassFile:
    :param pMassDcmpFile:
    :param pMassDcmpSpkesFile:
    :return:
    """

    if isinstance(pMassFile,str):
        if os.path.isfile(pMassFile):
            # Load in the data to see what the hell it is
            inputDataFile = open(pMassFile, "rb")
            dataOut = pickle.load(inputDataFile)
            inputDataFile.close()
            # Now we need to write a description of what the data is! [dataOut, xmin, xmax, vmin, vmax, amin, amax, dt, tmax]
            numbOfTrials = len(dataOut[0]) # This is the number of individual trails stored
            xmin = dataOut[1]
            xmax = dataOut[2]
            vmin = dataOut[3]
            vmax = dataOut[4]
            amin = dataOut[5]
            amax = dataOut[6]
            dt = dataOut[7]
            tmax = dataOut[8]
            descpFName = pMassFile[:-4]+'DESC.txt'
            file = open(descpFName,'w')
            file.write("1-D Moving Mass Data Description\n")
            file.write("\nData here is generated from the genAndSaveMoving1DMassData() method\n")
            file.write("\nFileName... "+pMassFile)
            file.write("\niterations. "+str(numbOfTrials))
            file.write("\nxmin....... "+str(xmin))
            file.write("\nxmax....... "+str(xmax))
            file.write("\nvmin........"+str(vmin))
            file.write("\nvmax........"+str(vmax))
            file.write("\namin........"+str(amin))
            file.write("\namax........"+str(amax))
            file.write("\ndt.........."+str(dt))
            file.write("\ntmax........"+str(tmax))
            file.close
    if isinstance(pMassDcmpFile, str):
        if os.path.isfile(pMassDcmpFile):
            # Load in teh data to see what it has it there
            inputDataFile = open(pMassDcmpFile, "rb")
            dataOut = pickle.load(inputDataFile)
            inputDataFile.close()
            # Now we need to write a description of what the data is! [segmentedTrials, gCenters, b, dataFile, xStart, xStop, nGaus, bInitial, dt]
            numbOfTrials = len(dataOut[0])
            originatingDataFile = dataOut[3]
            numbGaus = dataOut[6]
            xStart = dataOut[4]
            xStop = dataOut[5]
            GausVarb = dataOut[2]
            dt = dataOut[8]
            descpFName = pMassDcmpFile[:-4]+"DESC.txt"
            file = open(descpFName,'w')
            file.write("1-D Moving Mass Decomposed Data Description\n")
            file.write("\nData here is generated from the decomposeMoving1DMassData() method\n")
            file.write("\nFileName......."+pMassDcmpFile)
            file.write("\nInputDataFile.."+originatingDataFile)
            file.write("\niterations....."+str(numbOfTrials))
            file.write("\nNumbOfGauss...."+str(numbGaus))
            file.write("\nxStart........."+str(xStart))
            file.write("\nxStop.........."+str(xStop))
            file.write("\nGausVarbl_b...."+str(GausVarb))
            file.write("\ndt............."+str(dt))
            file.close()
    if isinstance(pMassDcmpSpkesFile, str):
        if os.path.isfile(pMassDcmpSpkesFile):
            # Load in teh data to see what it has it there
            inputDataFile = open(pMassDcmpSpkesFile, "rb")
            dataOut = pickle.load(inputDataFile)
            inputDataFile.close()
            # Now we need to write a description of what the data is! [segmentedSpikesList, dataFile, spikeGenType, analgSingnalScaling]
            numbOfTrials = len(dataOut[0])
            originatingDataFile = dataOut[1]
            spikeGenType = dataOut[2]
            analogScl = dataOut[3]
            descpFName = pMassDcmpSpkesFile[:-4]+"DESC.txt"
            file = open(descpFName,'w')
            file.write("1-D Moving mass Decomposed and Turned to Spikes Data Description\n")
            file.write("\nData here is generated from the decompedToSpikes1DMassData\n")
            file.write("\nFileName........."+pMassDcmpSpkesFile)
            file.write("\nInputDataFile...."+originatingDataFile)
            file.write("\nIterations......."+str(numbOfTrials))
            file.write("\nSpike-Gen-Type..."+spikeGenType)
            file.write("\nSignal-Scaling..."+str(analogScl))
            file.close()






########################################################################################################################
# RUNNING METHODS
########################################################################################################################
# Some the methods will be used to create acutal data and some will be used for testing. The differentiation will be
# given.

if __name__ == "__main__":


    ####################################################################################################################
    # TESTING METHODS
    ####################################################################################################################

    # print('Begin testing ... ... \n\n')

    # CREATE AND SAVE SOME 1D MOVING MASS DATA
    #genAndSaveMoving1DMassData()

    # PLOT THE CREATED AND SAVED DATA
    # We simply want to make sure that we are making what we are expecting.
    #loadAndPlot1DMassData()

    # DECOMPOSE 1D MOVING MASS DATA
    #decomposeMoving1DMassData()

    # PLOT THE 1D MOVING MASS DATA AND THE DECOMPOSITION OF SAID DATA
    #loadAndPlotDecomp1DMassData()

    # TURN DECOMPOSED 1D MOVING MASS DATA INTO SPIKES
    #decompedToSpikes1DMassData()

    # PLOT THE SPIKED VERSIONOF THE DECOMPOSED DATA FOR ACCURACY OF THE APPROACH
    #loadAndPlotDecomp1DMassData()

    # CREATE TXT FILE DESCRIPTIONS
    createDataDescriptionTxtFile(pMassFile='movingPointMassData/pointMassData000.pkl',pMassDcmpFile='movingPointMassData/pointMassDataDecmp000.pkl',pMassDcmpSpkesFile='movingPointMassData/pointMassDataDecmpSpikes000.pkl' )

    ####################################################################################################################
    # GENERATE EXPERIMENTAL DATA
    ####################################################################################################################
    # Here we actually will generate the data used in the experiments, we will need quite a bit of data. Note when
    # saving data the method zfill() seems pretty handy!

    # GENERATE 1D POINT MASS DATA
    #numberOfTrials = 10
    #topDirectory = 'movingPointMassData/'
    #subDirectory = 'dataSet_00/'
    #pMassDataName = 'pMassData'
    #pMassDcmpDataName = 'pMassDataDcmp'
    #pMassDcmpSpksDataName = 'pMassDcmpSpksData'
    #fileType = '.pkl'
    #First Step is to make the 1D Point Mass Data!
    #saveName = topDirectory+subDirectory+pMassDataName+fileType


