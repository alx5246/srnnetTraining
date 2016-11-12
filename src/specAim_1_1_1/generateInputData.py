# A.Lonsberry
# October 2016
#
# DESCRIPTION
#   Here I am going to have methods to create and save data for this experiment. The methods below are VERY hardcoded,
#   but rely on other methods that are generic. They are hardcoded so we can generate very specific data (original data,
#   decomposed original data, and then spiking from analog)

import matplotlib.pyplot as plt

# I had a difficult time getting things to import, so now I have however, go about fixing this by going up a directoy
# which I think I can do in the following way.
import os
import sys
import pickle
sys.path.append(os.path.abspath("../"))

# Call the python file that has the 1D mass simulator in it.
from dataGenerators import movingPointMass1D
from spikeDataGenerators import decompDimensions

########################################################################################################################
# METHODS FOR CREATING 1D MOVING MASS DATA (SMOOTH, DECOMPOSED, SPIKING)
########################################################################################################################

def genAndSaveMoving1DMassData(saveName='movingPointMassData/pointMassData000.pkl'):
    """
    DESCRIPTION
    Here we call the methods repeatedly to make instances of the 1D moving mass. Of these n examples, we then put them
    into a list, and then save them in the appropriate folder.
    :return: N/A (we save data to a file), the list saved is [listOfTrials, xmin, xmax, vmin, vmax, amin, amax, dt, tmax]
    """
    #How many iterations we want to include
    Iterations = 10
    xmin = 0.0
    xmax = 5.0
    vmin = 1.0
    vmax = 5.0
    amin = 1.0
    amax = 2.0
    dt   = .001
    tmax = 10.0
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
    trackingData = pickle.load(inputFile)
    inputFile.close()
    trackedArrays = trackingData[0]
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
    toSave = [segmentedTrials, gCenters, b, dataFile]
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
    # Fill in the

    # Load the data back (this is the decomposed version of the 1D moving mass data)
    inputDataFile = open(dataFile, "rb")
    dataOut = pickle.load(inputDataFile)
    inputDataFile.close()
    segmentedTrials = dataOut[0] # list of numpy.arrays(time-steps X numb. decomps) of the decomposed 1D signal
    gCenters = dataOut[1]  # The centers of the Gaussaians
    b = dataOut[2]
    origFileName = dataOut[3]

    # I need to load the original filename to get the value of 'dt'
    inputDataFile = open(origFileName, "rb")
    origDataOut = pickle.load(inputDataFile)
    inputDataFile.close()
    dt = origDataOut[7]


########################################################################################################################
# RUNNING METHODS
########################################################################################################################
# Some the methods will be used to create acutal data and some will be used for testing. The differentiation will be
# given.

if __name__ == "__main__":

    print('Begin testing ... ... \n\n')

    # CREATE AND SAVE SOME 1D MOVING MASS DATA
    #genAndSaveMoving1DMassData()

    # PLOT THE CREATED AND SAVED DATA
    # We simply want to make sure that we are making what we are expecting.
    #loadAndPlot1DMassData()

    # DECOMPOSE 1D MOVING MASS DATA
    #decomposeMoving1DMassData()

    # PLOT THE 1D MOVING MASS DATA AND THE DECOMPOSITION OF SAID DATA
    #loadAndPlotDecomp1DMassData()