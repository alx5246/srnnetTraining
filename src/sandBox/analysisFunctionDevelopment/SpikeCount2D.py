# S. Pickard
# June 2016
# (Using python 3.4.4, and brian2 2.0rc)
#
# DESCRIPTION
#   Function/method that produces an output vector of spike counts for each neuron in the network


import numpy as np

def spike_count2D(spikeTime, start, stop, dt):
    """
    FUNCTION DESCRIPTION
        This function takes as input a 2-D numpy.array of spike times, and outputs a spike count vector

    :param spikeTime: 2D numpy.array, units are seconds, neuron spike times stored in an numpy.array
    :param start: Start time of run in units of seconds
    :param stop: End time of run in units of seconds
    :param dt: user defined time subinteval over which to count spike
    """

    # Spike time turned into a numpy array
    spikeTime = np.array(spikeTime)
    # print ('spikeTime: ', spikeTime)

    # Creat interval array - intervals in which to break up the time array - sub time interval array
    duration = stop - start                             # Total run time
    n = duration / dt                                   # How many subintervals from time horizon results from user defined interval
    splitInterval = np.linspace(0, duration, n + 1)     # create numpy array of subinterval over which to count spikes
    # print ('splitInterval: ', splitInterval)

    #Setup counters for loop

    row=0       #Row iterator
    col=0       #iterator for columns
    j=0         #index for splitInterval array.
    k = 0       #index for new matrix that will store the grouped values from the split time array

    counter = 0
    SpikeCount = [[] for m in range(len(spikeTime))]

    for row in range(len(spikeTime)+1):
        if row < (len(spikeTime)):
            for col in range(len(splitInterval) + len(spikeTime[row])-1):
                if (col==0) and (spikeTime[row][col] == splitInterval[0]):
                    counter += 1
                    col += 1

                    #if/else statements to determine which interator gets incremented
                    if k < (len(spikeTime[row])):
                        k += 1
                    else:
                        j += 1

                elif (k < len(spikeTime[row])) and (spikeTime[row][k] > splitInterval[j]) and (spikeTime[row][k] <= splitInterval[j + 1]):
                    counter += 1
                    #if/else statements to determine which interator gets incremented
                    if k < (len(spikeTime[row])):
                        k += 1

                else:
                    SpikeCount[row].append(counter)
                    j += 1
                    if (k==len(spikeTime[row])):
                        #if j is less than the last time element
                        if (j<len(splitInterval)-1):
                            counter = 0
                            # if this is the last interval
                            if (j == len(splitInterval)-2):
                                j = 0
                                k = 0
                                counter = 0
                                SpikeCount[row].append(counter)
                                break
                        else: # j is the last time element
                            j = 0
                            k = 0
                            counter = 0
                    elif (j == (len(splitInterval)-1)):
                        j = 0
                        k = 0
                        counter = 0
                        break
                    else:
                        counter = 0

        else:
            if row < len(spikeTime):
                counter+=1
                SpikeCount[row].append(counter)
                counter = 0
                row += 1
                j = 0
                k = 0
            else:
                break

    return (SpikeCount, splitInterval)
