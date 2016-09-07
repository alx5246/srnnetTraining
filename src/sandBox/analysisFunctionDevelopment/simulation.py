# A. Lonsberry
# May 2016
# (Using python 3.4.4, and brian2 2.0rc)
#
# DESCRIPTION
#   Creating data, this will simply create data and save (using pickle) from a simple network. The data will be saved
#   in the folder 'savedData/'. The saved data (see code at bottom of the document) will be the arrays of data taken
#   from the simulation, specifically from brian2's 'StateMonitor' class.
#
# NETWORK DESCRIPTION
#   -Mostly taken from Carlson et al. "ramp" example (2013), some differences and addtions are made to equations
#   -100 input Poisson neurons
#   -1 RS Izhicevhich neuron
#   -STDP + Homeostatis (Carlson 2013) with some slight adaptations and alterations
#

from brian2 import *
import pickle
import numpy

import os
import shutil
import sys
import time as Time
import tarfile

currentDir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(currentDir)


class SimObject(object):
    def __init__(self, brianInternalDir=None):
        if brianInternalDir is not None:
            if os.path.exists(brianInternalDir):
                shutil.rmtree(brianInternalDir)
            os.makedirs(brianInternalDir)
        startTime = Time.clock()

        # this is different from the earlier simulation. brian2 will (by default) build the
        # network when a call to run() is made. brian2 is made so that only ONE network can
        # be built (rebuilding/making a new one is not allowed...Note that the term to refer
        # the the SINGLE instance of the network is "singleton"), so I changed the command
        # to NOT build the singleton when a call to run() is made (build_on_run=False).

        # What this means is that we need to make a call to build the singleton BEFORE
        # a call to run() is made. This is done at the end of the constructor.
        set_device('cpp_standalone', directory=os.path.join(currentDir, 'compiledCppSim'),
                   build_on_run=False)

        start_scope()  # Not sure why I put this in here...
        # With vlaues for (RS) type neuron used.
        a = .02  #
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
        # Setting up Poisson Neurons
        ########################################################################################################################

        P = PoissonGroup(100, numpy.arange(100) * .2 * Hz + .2 * Hz)

        ########################################################################################################################
        # STDP + Homeostatis Synapses (Carlson 2013)
        ########################################################################################################################

        # Synapse constants and parameters, mostly like those in (Carlson 2013)
        taupre = .020 * second
        taupost = .060 * second
        Apre = .0002
        Apost = -.000066
        R_tar = 6.6425  # This value is calculated for the low-pass filter and using 10Hz
        tauR = 1 * second

        # Synapse dynamics and models
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

        # Set the weights equivalently
        S.w = 2.0

        ########################################################################################################################
        # Monitors
        ########################################################################################################################

        # I will make the recording operate at a different rate by changing the time-step in the monitor (change "dt" value)

        # we need a SimObject instance to "remember" each monitor
        # so that in SimObject.execute(), the values can be pickled.
        # we accomplish this by making each monitor a "class field,"
        # (denoted by the self.<field_name>)
        self.SynMonitor = StateMonitor(S, ['w', 'R_hat', 'K'], record=numpy.arange(100), dt=.500 * ms)
        self.PoiMonitor = SpikeMonitor(P)
        self.NeuMonitor1 = SpikeMonitor(G)
        self.NeuMonitor2 = StateMonitor(G, ['v', 'u', 'L'], record=0, dt=.500 * ms)

        # this is the call to build the network. To do so, we need to get the device
        # that we set with the set_device() function. We can do this be calling the
        # get_device() function, because build() is actually a method (member function)
        # of whatever object brain2 is internally representing a network.
        get_device().build(directory=brianInternalDir)
        endTime = Time.clock()
        # print("Took %ss to construct Network using brian2" % (endTime - startTime))

    # a method to run the simulation on the network made in the constructor.
    def execute(self, alpha, beta, T, gamma, time, timestep, iteration, fileName):
        startTime = Time.clock()
        defaultclock.dt = timestep * usecond  # Set the simulation clock
        run(time * second)

        ########################################################################################################################
        ########################################################################################################################
        # SAVING
        ########################################################################################################################
        ########################################################################################################################
        # Now we want to save our results with the pickle function, IMPORTANTALY, brian has all these units floating around, we
        # need to get rid of them to use the pickle stuff, thus you will see things like "...numpy.asarray()..." to recast
        # and get rid of the units. (See brian2 online docs, http://brian2.readthedocs.io/en/2.0rc1/user/units.html)
        #
        # Also, by trial-and-error, I have found some ways to access the units of things....
        #   -If I have an array z for instance and I multiple by some unit like ms, I can get the unit by calling a.dimensions
        #   -If i have a statemonitor, called M, and I want to find the units of time for the inherent "t" array, I can call
        #    M.t.unit.
        # However this shit is weird and not easy to use (which things have which attributes). However brian has some built in
        # methods to check units by using the 'get_units(obj)' function.



        ########################################################################################################################
        # SAVING Version 1. (State-Monitor Output)
        ########################################################################################################################
        # Here I am saving each individual numpy array in an individual file. This works as long as I get rid of the brian units
        # that are associated with the array that is assigned by brian2.

        # Some print statements that I used to help debug the code that was not working.
        # print("The type for SynMonitor.t is: ",type(SynMonitor.t))
        # print("The type for NeuMonitor1.t is: ",type(NeuMonitor1.t))
        # print("The type for SynMonitor.t[0] is: ",type(SynMonitor.t[0]))
        # print("The type for SynMonitor.w is: ",type(SynMonitor.w))
        # print("The type for NeuMonitor1.t is: ",type(NeuMonitor1.t))
        # print("NeuMonitor1.t is \n", NeuMonitor1.t)
        # print("SynMonitor.t is \n", SynMonitor.t)
        # print("The type for numpy.asarry(NeuMonitor1.t) is: ", type(numpy.asarray(NeuMonitor1.t)))

        # print("\n....Saving data, version 1 ......")
        global currentDir
        saveDir = os.path.join(currentDir, fileName)
        os.makedirs(saveDir)

        with open(os.path.join(saveDir, "syn_w.pickl"), 'wb') as file:
            pickle.dump(numpy.asarray(self.SynMonitor.w), file)
        with open(os.path.join(saveDir, 'syn_w_units.pickl'), 'wb') as file:
            pickle.dump(get_unit(self.SynMonitor.w), file)

        with open(os.path.join(saveDir, "syn_t.pickl"), 'wb') as file:
            pickle.dump(numpy.asarray(self.SynMonitor.t), file)
        with open(os.path.join(saveDir, 'syn_t_units.pickl'), 'wb') as file:
            pickle.dump(get_unit(self.SynMonitor.t), file)

        with open(os.path.join(saveDir, "syn_r_hat.pickl"), 'wb') as file:
            pickle.dump(numpy.asarray(self.SynMonitor.R_hat), file)
        with open(os.path.join(saveDir, 'syn_r_hat_units.pickl'), 'wb') as file:
            pickle.dump(get_unit(self.SynMonitor.R_hat), file)

        with open(os.path.join(saveDir, "syn_k.pickl"), 'wb') as file:
            pickle.dump(numpy.asarray(self.SynMonitor.K), file)
        with open(os.path.join(saveDir, 'syn_k_units.pickl'), 'wb') as file:
            pickle.dump(get_unit(self.SynMonitor.K), file)

        with open(os.path.join(saveDir, "neu2_v.pickl"), 'wb') as file:
            pickle.dump(numpy.asarray(self.NeuMonitor2.v), file)
        with open(os.path.join(saveDir, 'neu2_v_units.pickl'), 'wb') as file:
            pickle.dump(get_unit(self.NeuMonitor2.v), file)

        with open(os.path.join(saveDir, "neu2_u.pickl"), 'wb') as file:
            pickle.dump(numpy.asarray(self.NeuMonitor2.u), file)
        with open(os.path.join(saveDir, 'neu2_u_units.pickl'), 'wb') as file:
            pickle.dump(get_unit(self.NeuMonitor2.u), file)

        with open(os.path.join(saveDir, "neu2_l.pickl"), 'wb') as file:
            pickle.dump(numpy.asarray(self.NeuMonitor2.L), file)
        with open(os.path.join(saveDir, 'neu2_l_units.pickl'), 'wb') as file:
            pickle.dump(get_unit(self.NeuMonitor2.L), file)

        with open(os.path.join(saveDir, "poi_t.pickl"), 'wb') as file:
            pickle.dump(numpy.asarray(self.PoiMonitor.t), file)
        with open(os.path.join(saveDir, 'poi_t_units.pickl'), 'wb') as file:
            pickle.dump(get_unit(self.PoiMonitor.t), file)

        with open(os.path.join(saveDir, "poi_i.pickl"), 'wb') as file:
            pickle.dump(numpy.asarray(self.PoiMonitor.i), file)
        with open(os.path.join(saveDir, 'poi_i_units.pickl'), 'wb') as file:
            pickle.dump(get_unit(self.PoiMonitor.i), file)

        with open(os.path.join(saveDir, "neu1_t.pickl"), 'wb') as file:
            pickle.dump(numpy.asarray(self.NeuMonitor1.t), file)
        with open(os.path.join(saveDir, 'neu1_t_units.pickl'), 'wb') as file:
            pickle.dump(get_unit(self.NeuMonitor1.t), file)

        with open(os.path.join(saveDir, "neu1_i.pickl"), 'wb') as file:
            pickle.dump(numpy.asarray(self.NeuMonitor1.i), file)
        with open(os.path.join(saveDir, 'neu1_i_units.pickl'), 'wb') as file:
            pickle.dump(get_unit(self.NeuMonitor1.i), file)

        with tarfile.open(os.path.join(currentDir, fileName + ".tar.gz"), "w:gz") as tarFile:
            tarFile.add(saveDir, arcname=fileName)
        shutil.rmtree(saveDir)
        # Now that I have saved the data I want to make sure I can also open it and it is what I expect!
        # print("Testing saved data, opening data just saved ... ")

        # For testing, we try opening the files using the pickle module/library
        # inputFile = open(os.path.join(currentDir, "savedData_0/netOutput0_PoiNeu_SpikesTimes.pkl"),"rb")
        # spikeTimes = pickle.load(inputFile)
        # spikeTimesUnits = pickle.load(inputFile)
        # inputFile.close()
        # print(spikeTimes)
        # print(spikeTimesUnits)



        # Create Dictionary
        endTime = Time.clock()
        # print("\tTook %ss to run simulation and record data" % (endTime - startTime))
        return os.path.join(currentDir, fileName + ".tar.gz")
