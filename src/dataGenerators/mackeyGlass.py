# A. Lonsberry
# September 2016
#
# DESCRIPTION
#   The Mackey-Glass equation is a nonlinear time delay DE that displays time-varying chaotic dynamics.
#
#   dx/dt = (beta*x(t-tau))/(1+x(t-tau)^n) - eta*x(t) f
#
#   An example in MATLAB uses beta = .2, eta=.1, tau=17, n=10, and x(0) = 1.2
#
#   This work here remains under construction.
#
# PROBLEMS
#   Some time the discretization is not sufficient, I need to be able to well handle this. This


import numpy
import matplotlib.pyplot as plt

def mackeyGlass(beta, eta, n, tau, x0, timeSteps):
    """
    DESCRIPTION
        Create data using mackey-glass time-delay differential equation

    PROBLEMS
        Some times I need a better stepping as the numerical integration can end up being too coarse, but when I
        simulate with a finer dt, then it fails to work... need to find a fix to this.

    :param beta:
    :param eta:
    :param n:
    :param tau:
    :param x0:
    :param dt:
    :param time:
    :return:
    """

    timeValues = numpy.arange(0., timeSteps, 1.) #Last step is NOT included
    print(numpy.shape(timeValues))
    output = numpy.ones(numpy.shape(timeValues)[0])*x0
    for i in range(tau,numpy.shape(output)[0]):
        output[i] = output[i-1] + ((beta*output[i-tau])/(1+output[i-tau]**n) - eta*output[i-1])

    return output


########################################################################################################################
# UNIT TESTING
########################################################################################################################

if __name__ == "__main__":

    plt.figure(1)
    mgValues = mackeyGlass(2, 1, 9.65, 2, 1.2, 100, 100)
    plt.plot(mgValues)

    plt.figure(2)
    plt.plot(mgValues[1:],mgValues[:-1])

    plt.show()
