# A. Lonsberry
# September 2016
#
# DESCRIPTION
#   I may want to simulate a moving object and its position vs. time. These methods are primarily used for generating
#   deterministic data for experimental testing.

import numpy
import random
import matplotlib.pyplot as plt

def randomMassMovement(xmin, xmax, v0min, v0max, a0min, a0max, dt, maxTime):
    """
    DESCRIPTION
    A single point-mass is moved in 1D given the input data.
    :param xmin:
    :param xmax:
    :param vmin: must be >= 0
    :param vmax: must be >= 0
    :param amin: must be >= 0
    :param amax: must be >= 0
    :param dt:
    :param maxTime: the maximum amount of time we record from
    :return: [positionArray, timeArray]
    """

    #Checking inputs
    if v0min >= 0. and v0max > v0min and a0min >= 0. and a0max > a0min and dt>0.0:

        # Randomly choose our initial condition (either edge of the domain bounds)
        if random.random() >.5:
            xCurrent = xmin
        else:
            xCurrent = xmax

        # Randomly choose initial value of velocity (make sure the initial velocity should push us through the domain)
        v0 = random.random()*(v0max-v0min) + v0min
        # Choose sign on velocity, based of course on the initial starting point of x0, as we want the point mass to travel
        # through the domain of interest.
        if xCurrent == xmax:
            v0 *= -1.

        # Now randomly choose value of acceleration and the sign of acceleration
        a0 = random.random()*(a0max-a0min) + a0min
        if random.random() > .5:
            a0 *= -1.

        # Instead of closed form, we fucking just simulate the particles position until it leaves the domain!
        # Create time-array, NOTE, because of rounding errors we have dumped numpy.arange() so not I am using
        # numpy.linspace, but have to be careful about including or excluding last points
        timeArray = numpy.linspace(0,maxTime,num=int(maxTime/dt+1))
        positionArray = numpy.zeros(int(maxTime/dt+1))
        positionArray[0] = xCurrent
        index = 1
        while xCurrent>=xmin and xCurrent<=xmax and index<=int(maxTime/dt):
            xCurrent = positionArray[0] + v0*timeArray[index] + .5*a0*(timeArray[index]**2.0)
            positionArray[index] = xCurrent
            index += 1

        #Truncate
        positionArray = positionArray[0:index]
        timeArray = timeArray[0:index]

        return [positionArray, timeArray]

    else:
        #Input constraints not met, return empty numpy array
        return [numpy.array([]), numpy.array([])]


########################################################################################################################
#                                              Internal Unit Testing
########################################################################################################################

if __name__ == "__main__":

    #Simply see if the function is working as expected to work.
    funcOuts = randomMassMovement(0., 5., 1.0, 5.0, 1.0, 2.0, .01, 10)
    plt.figure(1)
    plt.plot(funcOuts[1],funcOuts[0])
    plt.show()
