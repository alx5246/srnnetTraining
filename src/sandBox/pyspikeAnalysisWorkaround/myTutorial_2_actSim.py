# A. Lonsberry
# June 2016
# (Using python 3.4.4, and brian2 2.0rc, pyspike (uncompiled) )
#
# DESCRIPTION
#   A demonstration of using pyspike with brian2 simulation.
#
#   The brian2 simulation is similar to that of the ramp simulation by Carlson et al. (2013) where we have 100 Poisson
#   type neuron's firing into a singl LI&F neuron with STPD + Homeostatis learning rule.
#
#   For speedy computatoin, we also use brian2's standalone C++ compliaton option below.
#
#   After the simulation, we plot out some stuff using pyspike.
#
# REFERNCES
#   The pyspike module, was developed using rather newly published algorithms. To see a generic overview of these
#   approaches see the scholarpedia articles:
#   "Measures of spike train synchrony" : http://www.scholarpedia.org/article/Measures_of_spike_train_synchrony
#   "SPIKE - distance" : http://www.scholarpedia.org/article/SPIKE-distance
#
#   The Original papers wehere these measures have been develeoped are enumerate as follows:
#   (1) ISI-distance : Kreuz T, Haas JS, Morelli A, Abarbanel HDI, Politi A (2007a). "Measuring spike train synchrony". J Neurosci Methods 165:151–161.
#   (2) SPIKE - distance : Kreuz T, Chicharro D, Houghton C, Andrzejak RG, Mormann F (2013). "Monitoring spike train synchrony". JNeurophysiol 109:1457-72.
#   (3) Spike synchronization : ???
#

from brian2 import *
import numpy
import pyspike as spk

########################################################################################################################
# Setting up Simultation Comilation (brian2)
########################################################################################################################
# (brian2's 'standalone' compilation option is fast, thus we want to compile into c++ before running

# Nothing fancey where, we want to build on 'run'
set_device('cpp_standalone', directory='compiledCppSim')

########################################################################################################################
# Setting up General Simultation (brian2)
########################################################################################################################

start_scope()
defaultclock.dt = 100. * usecond

########################################################################################################################
# Setting up Izh. Neuron Model(s) (brian2)
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

# Make Neuron Group 1, where we have the single output Izh. type neuron of "RS" type
G = NeuronGroup(1, model=nurEqs, threshold='v > 30.', reset=resetEqs, refractory=2 * ms, method='euler')

########################################################################################################################
# Setting up Poisson Neurons (brian2)
########################################################################################################################

P = PoissonGroup(100, numpy.arange(100) * .2 * Hz + .2 * Hz)

########################################################################################################################
# STDP + Homeostatis Synapses (Carlson 2013) (brian2)
########################################################################################################################

# Synapse constants and parameters, mostly like those in (Carlson 2013)
taupre  = .020 * second
taupost = .060 * second
Apre    = .0002
Apost   = -.000066
R_tar   = 6.6425  # This value is calculated for the low-pass filter and using 10Hz
tauR    = 1 * second
beta    = 1
alpha   = .4
gamma   = 50
T       = 5

#Synapse dynamics and models
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

# Set the weights equivelently
S.w = 2.0

########################################################################################################################
# Monitors (brian2)
########################################################################################################################

# I will make the recording operate at a different rate by changing the time-step in the monitor (change "dt" value)
SynMonitor  = StateMonitor(S, ['w', 'R_hat', 'K'], record=numpy.arange(100), dt=.500*ms)
PoiMonitor  = SpikeMonitor(P)
NeuMonitor1 = SpikeMonitor(G)
NeuMonitor2 = StateMonitor(G, ['v', 'u', 'L'], record=0, dt=.500*ms)

########################################################################################################################
# RUN (brian2)
########################################################################################################################
run(1*second)

########################################################################################################################
# PuSpike ANALYSIS (pyspike)
########################################################################################################################
# Pyspike is a small module with limited functionality. Here we plot out the some of the results of processing with
# pyspike. In particular by using their ISI-distance, SPIKE-distance, and Spike-synchrony measures.

# Pull out the spike times (pull out the numpy arrays) from the spike monitors. Brian2's spike monitor class is unique
# to brian2, thus you cannot simpley give that object to any of pyspike's methods and expect it to handle the objects.
# What we do to go from SpikeMonitor to something pyspike can handle is firstly pull the numpy array's out of the Spike
# Monitors. Secondly then we take these arrrays and make "SpikeTrain" objects. SpikeTrain is a class in pyspike.
# Thridly you give these "SpikeTrain" objects to functions in pyspike for analysis.

#Grab all the spike times and neurons indices from the "SpikeMonitor" objects
poiSpikeTimes = numpy.asarray(PoiMonitor.t)
poiSpikeNeuro = numpy.asarray(PoiMonitor.i)
neuSpikeTimes = numpy.asarray(NeuMonitor1.t)

#Make a list where we will put the set of "SpikeTrain" objects
spikiTimesObjList = []

#Fill the list with spike trains for each neuron in the simulation
for i in numpy.unique(poiSpikeNeuro):
    inds = numpy.argwhere(poiSpikeNeuro==i) #FInd tthe spike-times that correspond to neuron i
    inds = inds[:,0]                        #the variable "inds" is a 2D array, I just want a 1D array to give to the "SpikeTrain" class.
    spikiTimesObjList.append(spk.SpikeTrain(poiSpikeTimes[inds], edges=(0., 1.), is_sorted=True)) #Append the spikiTimesObjList with a "SpikeTrain" object

#Make a "spike_profile", which essentially is a simultaneous comparison of all the spike trains together simutaneously
spike_profile = spk.spike_profile(spikiTimesObjList)

#Make a "isi_profile", which essentially is a simultaneous comparison of all the spike trains together simutaneously
isi_profile = spk.isi_profile(spikiTimesObjList)

#Make a "spike_sync_profile", which essentially is a simultaneous comparison of all the spike trains together simutaneously
sync_profile = spk.spike_sync_profile(spikiTimesObjList)

#Now I want to pring out want these profiles that have just been created look like.
plt.figure(1)
x, y = spike_profile.get_plottable_data()
plt.plot(x, y, '--k')
plt.figure(2)
x, y = isi_profile.get_plottable_data()
plt.plot(x,y, '--k')
plt.figure(3)
x, y = sync_profile.get_plottable_data()
plt.plot(x,y, '--k')
plt.show()

