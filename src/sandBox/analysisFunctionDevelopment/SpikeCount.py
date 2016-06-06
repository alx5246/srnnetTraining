# S. Pickard
# June 2016
# (Using python 3.4.4, and brian2 2.0rc)
#
# DESCRIPTION
#   Function/method that produces an output vector of averaged spike counts across neurons per time step
#   This analysis was inspired by Lazar et al.(2009) Figure 5.b
#   Goal is to find how many spike occur in each sub time interval

import numpy as np

def spike_count(time, duration, dt):
    """
    FUNCTION DESCRIPTION
        This function takes as input a 1-D numpy.array of spike times, and outputs a spike count vector; the spike
        counts are averaged over a user defined interval

        This is an example function, and thus has some limitations. For example it can only handle 1D input. A more
        complex version of this function would likely be able to handle N-D input (spike trains from multiple neurons)
        simultaneously and would average the number of spikes over a user definesd interval and then average over the
        number of neurons so that each averaging interval, a scalar mean would be outputed. This function is not
        concerned with identifying how many spikes per interval are occurring for each neuron, but rather an averaged
        count across all neurons (i.e. the network) per interval

    :param time: 1D numpy.array, units are seconds, neuron spike times stored in an numpy.array
    :param duration: the time duration (in seconds) of the data collection run
    :param dt: user defined time subinteval over which to count spike
    """


    #Spike time turned into a numpy array
    time = np.array(time)

    #Interval array - intervals in which to break up the time array - sub time interval array
    n = duration/dt                                 #How many subintervals from time horizon results from user defined interval
    splitInterval = np.linspace(0, duration, n+1)   #create numpy array of subinterval over which to count spikes
    print splitInterval

    i=0     #inex for time array
    j=0     #index for splitInterval array.
    k=0     #index for new matrix that will store the grouped values from the split time array
    counter = 0     #counter variable to keep track of spike count for each subinterval through loop
    SpikeCount = []     #Initialize array to collect the number of spikes occuring wihtin each subinterval

    for i in xrange(len(time)):
        if (time[i] >= splitInterval[j]) & (time[i] <= splitInterval[j + 1]):
            counter += 1
            print ('counter: ', counter)
            # print ('time element: ', time[i])
            # print splitInterval[j]
            # print splitInterval[j + 1]
            i += 1
            print ('i: ', i)
        else:
            SpikeCount.append(counter)
            print ('Spike count: ', SpikeCount)
            # print "else statement"
            counter = 0
            counter += 1
            # k += 1
            i += 1
            j += 1
    SpikeCount.append(counter)
    print ('Spike count: ', SpikeCount)

    return SpikeCount
