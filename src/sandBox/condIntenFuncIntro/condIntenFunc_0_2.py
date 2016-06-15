# A. Lonsberry
# June 2017
#
# DESCRIPTION
#   This python file will have methods to run a brian2 neural simulation and outputs

from brian2 import *
import os
import numpy

def singleLayerSimulationNew(inputNetStruc, mainNetStruc, inputToMainConnections, inputToMainConnWeights, mainToMainConnections, mainToMainConnWeights, inputSpikes):
    '''
    DESCRIPTION
    This runs a new simulation with a Izh. RS neurons (slgihtly modified) along with a version of STDP + Homeostatsis
    developed in Carlson et al. 2013

    :param inputNetStruc:
    :param mainNetStruc:
    :param inputToMainConnections:
    :param mainToMainConnections:
    :param inputSpikes:
    :return:
    '''

    ####################################################################################################################
    # Setting up Simulation C++ Standalone Compilation (as per brain 2.0rc)
    ####################################################################################################################

    # Setting compilation for simply compiling inside a for-loop (this gets called inside a for-loop outside), and not
    # with anything multiprocessing.
    device.reinit()
    buildDirectory = 'standalone{}'.format(os.getpid())
    set_device('cpp_standalone', build_on_run=False)

    ####################################################################################################################
    # Setting up General Simulation
    ####################################################################################################################

    start_scope()
    defaultclock.dt = 100. * usecond

    ####################################################################################################################
    # Setting up Izh. Neuron Model(s)
    ####################################################################################################################

    # The Neuron model is nearly the same an Izh. (RS) neuron. The notable difference is that synapses not have
    # exponential current injectsions whcih is more typical of regular LI&F neuron models. This current will be
    # represented as a state-variable "L"

    # With vlaues for (RS) type neuron used.
    a = .02
    b = .2
    c = -65.0
    d = 8.0

    # I need time-constant values
    tau1 = 1000 * ms     #
    tau2 = 250 * ms      # This is the time-constant associated with current input

    # Neuron dynamics
    nurEqs = """
    dv/dt = (0.04 * (v**2) + 5.0*v + 140.0 - u + I + 5.*L)/tau1 : 1 (unless refractory)
    du/dt = (a * ( b * v - u ) ) / tau1 : 1
    dL/dt = -L/tau2 : 1
    I = 50.*(randn()+1) : 1
    """
    resetEqs = """
    v = c
    u = u + d
    """

    ####################################################################################################################
    # Setting up Main Network Bodya of neuonrs
    ####################################################################################################################

    # Neuron Group, where we have one neuron now for each different iteration of the values we care about. The number
    # of neurons in this group is define by the "mainNetStruc" array.

    G = NeuronGroup(mainNetStruc.shape[0], model=nurEqs, threshold='v > 30.', reset=resetEqs, refractory=2*ms, method='euler')

    ####################################################################################################################
    # Setting up Input Neuron Groupd
    ####################################################################################################################

    # We are giving the input neuron spike times as input into this method/fuction.

    inds  = inputSpikes[:, 0]
    times = inputSpikes[:, 1]
    times = times * second
    #print(inds.astype(int))
    #print(numpy.max(inds.astype(int)))
    P = SpikeGeneratorGroup(numpy.max(inds.astype(int))+1, inds.astype(int), times)

    ####################################################################################################################
    # STDP + Homeostatis Synapses
    ####################################################################################################################

    # The STDP + Homeostatis model here is largely a replication from Carlson et al. (2013).

    # Synapse dynamics are controlled by a set of parameters
    taupre  = .020 * second
    taupost = .060 * second
    Apre    = .0002
    Apost   = -.000066
    R_tar   = 6.6425  # This value is calculated for the low-pass filter and using 10Hz ... how did I calc. this? -AJL
    tauR    = 1. * second
    beta    = 0.05
    alpha   = 0.20
    gamma = 50.
    T = 5.

    # Make the synapes between input to main neuron group
    S = Synapses(P, G, model="""
                             dw/dt = (alpha * w * (1 - R_hat / R_tar) * K) / (1*second)  : 1 (clock-driven)
                             dapre/dt  = -apre / taupre : 1 (event-driven)
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

    # Make the synapses between main body to main body
    S2 = Synapses(G, G, model="""
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


    ####################################################################################################################
    # Connect synapses
    ####################################################################################################################

    # Connect all neurons in group 'P' to group 'G'
    S.connect(i = inputToMainConnections[:, 0], j = inputToMainConnections[:, 1])
    print(inputToMainConnections)
    if inputToMainConnWeights.size != 0:
        S.w = inputToMainConnWeights
    else:
        S.w = 10.0 #Initialize weights of the synapses



    # Connect all neurons in group 'G' to group 'G'
    S2.connect(i = mainToMainConnections[:, 0], j = mainToMainConnections[:, 1])
    if mainToMainConnWeights.size != 0:
        S2.w = mainToMainConnWeights
    else:
        S2.w = .95 #Initialize weights fo the synapses.

    ####################################################################################################################
    # Monitors
    ####################################################################################################################

    # I record all the spikes, also I want to record riring rates as well!

    M0 = SpikeMonitor(P)
    M1 = SpikeMonitor(G)

    M2 = PopulationRateMonitor(G)

    M3 = StateMonitor(G, ['I', 'v'], record=[0], dt=1*ms)

    M4 = StateMonitor(S, ['w'], record=[0], dt=100*ms)

    M5 = StateMonitor(S2, ['w'], record=[0], dt=100*ms)

    ####################################################################################################################
    # RUN
    ####################################################################################################################

    # Setting compilation for simpley compiling inside a for-loop (this gets called inside a for-loop outside), and not
    # with anything multiprocessing.

    runTime = numpy.max(times)
    run(runTime + 1.*second, report='text')
    device.build(directory=buildDirectory, compile=True, run=True, debug=False)

    ####################################################################################################################
    # RETURN
    ####################################################################################################################

    # Other than spike monitors, I would like to return the synaptic weights as well!
    outputList = [M0, M1, M2, S, S2, M3, M4, M5]

    return outputList

