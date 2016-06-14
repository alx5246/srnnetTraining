# S. Pickard
# June 2016
# (Using python 3.4.4, and brian2 2.0rc)
#
# DESCRIPTION
#   Function/method that produces an output vector of averaged spike counts across neurons per time step
#   This analysis was inspired by Lazar et al.(2009) Figure 5.b
#   Goal is to find how many spike occur in each sub time interval

import numpy as np

def spike_count(spikeTime, start, stop, dt):
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

    :param spikeTime: 1D numpy.array, units are seconds, neuron spike times stored in an numpy.array
    :param start: Start time of run in units of seconds
    :param stop: End time of run in units of seconds
    :param dt: user defined time subinteval over which to count spike
    """


    #Spike time turned into a numpy array
    spikeTime = np.array(spikeTime)
    print('Spike Times: ', spikeTime)

    #Creat interval array - intervals in which to break up the time array - sub time interval array
    duration = stop-start                           #Total run time
    n = duration/dt                                 #How many subintervals from time horizon results from user defined interval
    splitInterval = np.linspace(0, duration, n+1)   #create numpy array of subinterval over which to count spikes
    print ('split interval: ', splitInterval)

    ##Find length over which to iterate in for loop
    length_splitInt = len(splitInterval)
    print('length splitInterval: ', length_splitInt)
    length_time = len(spikeTime)
    print('length time: ', length_time)
    length = length_splitInt + ((length_time) - 2)
    print('length :', length)

    i=0                 #inex for time array
    j=0                 #index for splitInterval array.
    k=0                 #index for new matrix that will store the grouped values from the split time array
    counter = 0         #counter variable to keep track of spike count for each subinterval through loop
    SpikeCount = []     #Initialize array to collect the number of spikes occuring wihtin each subinterval

    for i in range(length):
        if (i == 0) and (spikeTime[0] == splitInterval[0]):
            counter += 1
            i += 1

            # Spot check
            print('if counter: ', counter)
            print('time element: ', spikeTime[k])
            print('splitInt: ', splitInterval[j], splitInterval[j + 1])
            print('i: ', i)
            print('if k: ', k)

            if k < (len(spikeTime) - 1):
                k += 1

                # Spot check
                print('iff k: ', k)
                print('iff counter: ', counter)
            else:
                j += 1

                # Spot check
                print('iff counter: ', counter)
                print(SpikeCount)
                print('iff j: ', j)

        elif (spikeTime[k] > splitInterval[j]) and (spikeTime[k] <= splitInterval[j + 1]):
            counter += 1
            i += 1

            # Spot check
            print('if counter: ', counter)
            print('time element: ', spikeTime[k])
            print('splitInt: ', splitInterval[j], splitInterval[j + 1])
            print('i: ', i)
            print('if k: ', k)

            if k < (len(spikeTime) - 1):
                k += 1

                # Spot check
                print('iff k: ', k)
                print('iff counter: ', counter)

            else:
                j += 1
                # Spot check
                SpikeCount.append(counter)
                print('iff counter: ', counter)
                print(SpikeCount)
                print('iff j: ', j)



        else:
            SpikeCount.append(counter)
            counter = 0
            j += 1
            i += 1

            # Spot Check
            print('else counter: ', counter)
            print(SpikeCount)
            print('time element: ', spikeTime[k])
            # print('splitInt: ', splitInterval[j], splitInterval[j + 1])
            print('else j: ', j)
            print('else i: ', i)
            print('else k: ', k)

    return SpikeCount
