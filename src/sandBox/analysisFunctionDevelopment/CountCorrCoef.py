# S. Pickard
# July 2016
# (Using python 3.4.4, and brian2 2.0rc)
#
# DESCRIPTION
#   Function/method that produces a scalar output of averaged correlation coefficents (of spike counts across neurons per
#   time subinterval) for a time duration
#   This analysis was inspired by Lazar et al.(2009) Figure 5.c
#   Goal is to find the average correlation between neurons' firing counts over a time duration (not subinterval)

import numpy as np

def count_Ave_CorrCoef(NeuCountArray):
    """
    FUNCTION DESCRIPTION
        This function takes as input a N-D numpy.array of spike times, and outputs a spike count vector; the spike
        counts are averaged over a user defined interval

    :param NeuCountArray: N-D numpy.array, units are seconds, neuron spike times for each neuron stored in an numpy.array
    """


    #Spike times array turned into a numpy array
    NeuCountArray = np.array(NeuCountArray)
    print(NeuCountArray)

    #Generate array of normalized correlation coefficients (makes a symmetrical matrix)
    CountCorrCoef = np.corrcoef(NeuCountArray, rowvar = True)
    print(CountCorrCoef)

    #Keep the upper triangle of symmetrical corr. coef. matrix and the diagnol of ones (the redunate information is converted to zeros)
    UpTriCorrCoef = np.triu(CountCorrCoef)
    print (UpTriCorrCoef)

    #Eliminate the zeros (redunant info from symmetrical matrix) from the symmetrical corr. coef matrix
    NoZero = np.extract(abs(UpTriCorrCoef) > 0, UpTriCorrCoef)
    print(NoZero)

    #Eliminate the ones from the symmetrical corr. coef matrix (eliminate the variances and keep covariances)
    CovariancesOnly = np.extract(abs(NoZero) < 1, NoZero)
    print(CovariancesOnly)

    #Average the covariances
    AveCov = np.mean(abs(CovariancesOnly))

    return (AveCov)