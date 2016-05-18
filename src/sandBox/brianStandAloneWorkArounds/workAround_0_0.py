#
# DESCRIPTION
#   I am showing a simple example of how you can use one simulation, to test a number of parameters. This essentially
#   boils down to running multiple networks simultaneoulsy.
#
#   In this example here I am running a bunch of Poisson inputs to a number of different neurons. The differences come
#   from the fact that synapse dynamics are different based on which target neuron the synapse connects into. Thus if I
#   have ten target neurons, I essentially make 10 different types of synapses, and get to run 10 different simulations
#   at the same time.

from brian2 import *
import matplotlib.pyplot as plt
import numpy
import numpy.matlib

#Setting up the C++ compilation, left as default which is build on run.
#set_device('cpp_standalone', directory='workAround_0_0_Comp')
set_device('cpp_standalone', directory='myTurorial_2_Izhecevich_2_Comp', build_on_run=False)

########################################################################################################################
# Setting up Multiple Simulations
########################################################################################################################
# What are we iterating over? This will effect how the simulation is generated
iterationParameter = "alpha"
iterationNumber = 10
iterationValues = numpy.linspace(start=0., stop=1., num=iterationNumber)

#reformat Values
iterationValues = numpy.matlib.repmat(iterationValues,1,100)
iterationValues = iterationValues[0,:]
#iterationValues = numpy.sort(iterationValues)

########################################################################################################################
# Setting up General Simultation
########################################################################################################################
start_scope()
defaultclock.dt = 100. * usecond

########################################################################################################################
# Setting up Izh. Neuron Model(s)
########################################################################################################################
#With vlaues for (RS) type neuron used.
a = .02
b = .2
c = -65.0
d = 8.0

#I need time-constant values
tau1 = 1000*ms
tau2 = 250*ms


#Neuron dynamics
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

#Make Neuron Group, where we have one neuron now for each different iteration of the values we care about
G = NeuronGroup(iterationNumber, model=nurEqs, threshold='v > 30.', reset=resetEqs, refractory=2*ms)

P = PoissonGroup(100, numpy.arange(100)*.2*Hz + .2*Hz)

########################################################################################################################
# STDP + Homeostatis Synapses
########################################################################################################################
#Synapse dynamics and group
taupre  = .020 * second
taupost = .060 * second
Apre    = .0002
Apost   = -.000066
R_tar   = 6.6425 #This value is calculated for the low-pass filter and using 10Hz
tauR = 1 * second
beta    = 1
#alpha   = 1
gamma   = 50
T       = 5
S = Synapses(P, G, model="""
                         dw/dt = (alpha * w * (1 - R_hat / R_tar) * K) / (1*second)  : 1
                         dapre/dt = -apre / taupre : 1 (event-driven)
                         dapost/dt = -apost / taupost : 1 (event-driven)
                         dR_hat/dt = -R_hat/tauR : 1
                         K = R_hat / ( T * abs( 1 - R_hat / R_tar) * gamma) : 1
                         alpha : 1
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
                         """)

#Connect all neurons in group 'P' to group 'G'
S.connect()

#Set the alpha values
S.alpha = iterationValues

#print(S.alpha[0,0])
#print(S.alpha[0,1])
#print(S.alpha[0,2])
#print(S.alpha[1,0])
#print(S.alpha[2,0])
#print(S.alpha[3,0])
#print(S.alpha)

#Set the weights
#vectWeights = numpy.linspace(start=1.0,stop=3.0,num=100)
#S.w = vectWeights #USE STRINGS IN COMPILED VERSION
S.w = 2.0

########################################################################################################################
# Monitors
########################################################################################################################

#M0 = SpikeMonitor(P)
M1 = StateMonitor(S,['R_hat'],record=[0,1,2,3,4,5,6,7])
M2 = StateMonitor(S,['w'],record=[0,1,2,3,4,5,6,7])
#M2 = StateMonitor(S,['alpha'],record=[0,100,200,300])
#M2 = StateMonitor(S,'w',record=numpy.arange(100))
#M2 = StateMonitor(S,['w','K',],numpy.arange(100), dt=500*ms)
#M3 = SpikeMonitor(G)
#M4 = PopulationRateMonitor(G)

########################################################################################################################
# RUN
########################################################################################################################

#run(500*second)
run(600*second, report='text')


########################################################################################################################
# Plot
########################################################################################################################

plt.subplot(611)
plt.plot(M1.t, M1.R_hat[0], '-b')
plt.subplot(612)
plt.plot(M1.t, M1.R_hat[1], '-b')
plt.subplot(613)
plt.plot(M1.t, M1.R_hat[2], '-b')
plt.subplot(614)
plt.plot(M1.t, M1.R_hat[3], '-b')
plt.subplot(615)
plt.plot(M1.t, M1.R_hat[4], '-b')
plt.subplot(616)
plt.plot(M1.t, M1.R_hat[5], '-b')

plt.figure(2)
plt.plot(M2.t, M2.w[0], '-b')
plt.subplot(612)
plt.plot(M2.t, M2.w[1], '-b')
plt.subplot(613)
plt.plot(M2.t, M2.w[2], '-b')
plt.subplot(614)
plt.plot(M2.t, M2.w[3], '-b')
plt.subplot(615)
plt.plot(M2.t, M2.w[4], '-b')
plt.subplot(616)
plt.plot(M2.t, M2.w[5], '-b')
plt.show()