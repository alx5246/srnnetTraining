#
# DESCRIPTION
#   Two parts: (orkAround_1_0.py and workAround_1_1.py.
#
#   In the other file we encapsulate the simulation within a function called 'runSimulation', that takes input
#   'alpha'. IMPORTANTLY at the top of the function we wet the device and all that so that the code gets compiled using
#   c++.
#
#   In the file here, workAround_1_1.py, we call this function over and over again inside a for-loop. Each time we make
#   a call to this function, there is overhead because we have to recompile each time. The overhead is worth it however
#   if the simulation has lots of dynamics.
#
#   An important fact is that in order for this encapsulation inside a function to work correctly you have to add a line
#   with 'device.reinit()'. I found this on the brian forum. The reason being is we want the compiler to forget about
#   prior iterations or prior calls of this simulation. Without the 'device.reinit()' this function will work once but
#   will fail the next time it is called.

import matplotlib.pyplot as plt
from workAround_1_0 import runSimulation
import numpy
from joblib import Parallel, delayed

########################################################################################################################
# Regular non-parallel version (Works fine)
########################################################################################################################
#
#alphaValues = numpy.array([0.,.2,.4,6,.8,1])
#
#
#for idx, alpha in enumerate(alphaValues):
#    M1 = runSimulation(alpha)
#    plt.subplot(alphaValues.shape[0],1,idx+1)
#    plt.plot(M1.t, M1.R_hat[0], '-b')
#
#plt.show()
#
########################################################################################################################
# Now try parallel version
########################################################################################################################
#
#Must use 'main' with windows, or this will span an infinite number of threads of something
if __name__ == '__main__':
    a = Parallel(n_jobs=4)(delayed(runSimulation)(i) for i in numpy.linspace(start=.2, stop=.8, num=4))
    print(a)
    print(type(a))