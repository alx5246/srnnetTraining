A. Lonsberry
May 2016

(a)nalysis (f)unction (d)evelopment README (afdREADME)

DESCRIPTION

    Here we will have an area dedicated to developing and testing analysis functions, for example in the case of examining
    firing-rates, firing-rate distributions, or spike train correlations etc.

FILES AND FOLDERS

    (1) createData_0_0.py and savedData_0/
    Short Description: creates data for use in testing analysis functions
    Long Description: All the files in this section of the sandBox/ folder are developed for the analysis of networks
    and their output. In order to verify these functions and that they are working correctly, we need to create data.
    The file, createData_0_0.py is a python that runs a simple simulation, mimicking that of the ramp test (Carlson
    2013), and then saves the data using pickle (pickle is a python module for saving and opening data) to the folder
    savedData_0/

    (2) firingRateDistrib_0_0.py and firingRateDistrib_0_0.py
    Short Description: example function/method to calculate a smooth firing-rate approximation given some spike-train.
    Long Description: One potentially important tool for visualization and analysis, is estimating the firing-rate
    of some neuron. The file firingRateDistrib_0_0.py, has the function to actually do the calculations. The file
    firingRateDistrib_0_1.py has code to call that function, and plot out the results so we can see the function is
    working as expected.

    (3) SpikeCount.py
    Short Description: spike count over time subintervals
    Long Description

    (4) CountPCA
    Short Description: principal compenent analysis script of spike times
    Long Description

    (5) spikeMonToMatrix.py
    Short Description: Find element positions for each neuron with the neuron index array (spikemon.i) generated from
    the spikmon class in the simulation. Once the elements for each neuron are found, these elements are then used to
    find the corresponding spike time for said neuron. The ouput is a series of arrays where each row is a neuron index
    and the columns are spike times

    (6) SpikeCount2D.py
    Short Description: spike count over time subintervals for a 2D vector input (spike times for each neuron)

    (7) CountCorrCoef.py
    Short Description: CountCorrCoef: takes in N-D array of neuron spike counts for each user defined subinterval (dt)
    over a total time duration (from script SpikeCount2D). Output is an averaged correlation coefficient for the time
    duration
