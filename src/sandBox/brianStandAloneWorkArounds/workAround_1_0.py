# A. Lonsberry
# May 2016
# (Using python 3.4.4, and brian2 2.0rc)
#
# DESCRIPTION
#   Two parts: this files (workAround_1_0.py) and workAround_1_1.py.
#
#   Here in the first file we encapsulate the simulation within a function called 'runSimulation', that takes input
#   'alpha'. IMPORTANTLY at the top of the function we wet the device and all that so that the code gets compiled using
#   c++.
#
#   In the second file workAround_1_1.py, we call this function over and over again inside a for-loop. Each time we make
#   a call to this function, there is overhead because we have to recompile each time. The overhead is worth it however
#   if the simulation has lots of dynamics.
#
#   An important fact is that in order for this encapsulation inside a function to work correctly you have to add a line
#   with 'device.reinit()'. I found this on the brian forum. The reason being is we want the compiler to forget about
#   prior iterations or prior calls of this simulation. Without the 'device.reinit()' this function will work once but
#   will fail the next time it is called.

from brian2 import *
import numpy
import numpy.matlib
import os

def runSimulation(alpha):

    ########################################################################################################################
    # Setting up Simultation Comilation
    ########################################################################################################################

    # Setting compilation for simpley compiling inside a for-loop (this gets called inside a for-loop outside), and not
    # with anything multiprocessing.
    #device.reinit()
    #buildDirectory = 'standalone{}'.format(os.getpid())
    #set_device('cpp_standalone', build_on_run=False)

    # Setting compilation for simpley compiling inside a for-loop while the for-loop is also includes parallel processing
    device.reinit()
    buildDirectory = 'standalone{}'.format(os.getpid())
    set_device('cpp_standalone', build_on_run=False)


    ########################################################################################################################
    # Setting up General Simultation
    ########################################################################################################################
    start_scope()
    defaultclock.dt = 100. * usecond

    ########################################################################################################################
    # Setting up Izh. Neuron Model(s)
    ########################################################################################################################
    # With vlaues for (RS) type neuron used.
    a = .02
    b = .2
    c = -65.0
    d = 8.0

    # I need time-constant values
    tau1 = 1000 * ms
    tau2 = 250 * ms

    # Neuron dynamics
    nurEqs = """
    dv/dt = (0.04 * (v**2) + 5*v + 140 - u + I + 5*L)/tau1 : 1 (unless refractory)
    du/dt = (a * ( b * v - u ) ) / tau1 : 1
    dL/dt = -L/tau2 : 1
    I : 1
    """
    resetEqs = """
    v = c
    u = u + d
    """

    # Make Neuron Group, where we have one neuron now for each different iteration of the values we care about
    G = NeuronGroup(1, model=nurEqs, threshold='v > 30.', reset=resetEqs, refractory=2 * ms, method='euler')

    P = PoissonGroup(100, numpy.arange(100) * .2 * Hz + .2 * Hz)

    ########################################################################################################################
    # STDP + Homeostatis Synapses
    ########################################################################################################################
    # Synapse dynamics and group
    taupre = .020 * second
    taupost = .060 * second
    Apre = .0002
    Apost = -.000066
    R_tar = 6.6425  # This value is calculated for the low-pass filter and using 10Hz
    tauR = 1 * second
    beta = 1
    # alpha   = 1 #Imported variable
    gamma = 50
    T = 5
    S = Synapses(P, G, model="""
                             dw/dt = (alpha * w * (1 - R_hat / R_tar) * K) / (1*second)  : 1 (clock-driven)
                             dapre/dt = -apre / taupre : 1 (event-driven)
                             dapost/dt = -apost / taupost : 1 (event-driven)
                             dR_hat/dt = -R_hat/tauR : 1 (clock-driven)
                             K = R_hat / ( T * abs( 1 - R_hat / R_tar) * gamma) : 1
                             """,
                 on_pre="""
                             L += w
                             apre += Apre
                             w += beta * K * apost
                             """,
                 on_post="""
                             apost += Apost
                             w += beta * K * apre
                             R_hat += 1
                             """,
                 method='euler')

    # Connect all neurons in group 'P' to group 'G'
    S.connect()

    # print(S.alpha[0,0])
    # print(S.alpha[0,1])
    # print(S.alpha[0,2])
    # print(S.alpha[1,0])
    # print(S.alpha[2,0])
    # print(S.alpha[3,0])
    # print(S.alpha)

    # Set the weights
    # vectWeights = numpy.linspace(start=1.0,stop=3.0,num=100)
    # S.w = vectWeights #USE STRINGS IN COMPILED VERSION
    S.w = 2.0

    ########################################################################################################################
    # Monitors
    ########################################################################################################################

    # M0 = SpikeMonitor(P)
    M1 = StateMonitor(S, ['R_hat'], record=[0], dt=1*ms)

    ########################################################################################################################
    # RUN
    ########################################################################################################################

    # Setting compilation for simpley compiling inside a for-loop (this gets called inside a for-loop outside), and not
    # with anything multiprocessing.
    #run(600 * second, report='text')
    #device.build(directory=buildDirectory, compile=True, run=True, debug=False)
    #return M1  # This works fine without multiprocess

    # Setting compilation for simpley compiling inside a for-loop while the for-loop is also includes parallel processing
    run(600 * second, report='text')
    device.build(directory=buildDirectory, compile=True, run=True, debug=False)
    #
    # If I give a simple dummy output this seems to work, but what about stuff I actually want out?
    #aRa = [1, 3, 4]
    #return aRa
    #
    # This does not work....
    #return M1 #apparently M1 is not picklable
    #print(type(M1.R_hat))
    #print(M1.R_hat)
    return M1.R_hat

