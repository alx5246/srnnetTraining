import datetime
import os
import sys


currentDir = os.path.dirname(os.path.realpath(__file__)) # returns a string
layerStripping = [".." for x in range(0, 4)] # a list of ".." strings ["..", "..", "..", ".."]

# RELATIVE TO THE DIRECTORY controller.py is in...currentDir is the string form
# of the file path to the directory that controller.py is in...modify it
dbPath = os.path.join(currentDir, *layerStripping) # equivalent of os.path.join(currentDir, "..", "..", "..", "..")

# let python know where to search
sys.path.extend([currentDir, dbPath])

import db
import simulation
import time as Time

startClockTime = Time.clock()

manager = db.DBManager(databaseName="carlsonOneNeuron")
manager.openCollection("STDParameterRange")

# I changed the time to 1000 seconds for you. It should still be running
# when you get back (estimated time is ~66000 seconds = 18.3 hours to run and this is a
# low ball estimate because I didn't account for setting up the network).
time = 1000   #seocnds
timestep = 100  #In units of microseconds
# gamma = 50
# beta = 1
# alpha = .4
# T = 5 # units of seconds
#Setup paramter lists
param1 = {
    "alpha": [float(x) / 10.0 for x in range(0, 11)],
}

param2 = {
    "beta": [float(x) / 10.0 for x in range(0, 110)],
}
param3 = {
    "gamma": [float(x) / 1.0 for x in range(0, 110)],
}
param4 = {
    "T": [float(x) / 1.0 for x in range(0, 110)],
}

print("MONGODB_URI: %s" % os.environ.get("MONGODB_URI"))

startTime = None
endTime = None
results = None

# Make sure to construct the network!! Only make ONE of these objects
sim = simulation.SimObject()


for alpha in param1["alpha"]:
    for beta in param2["beta"]:
        for gamma in param3["gamma"]:
            for T in param4["T"]:
                for simRep in range(10):
                    print("RUNNING SIMULATION WITH alpha: %s" % alpha)
                    print("RUNNING SIMULATION WITH beta: %s" % beta)
                    print("RUNNING SIMULATION WITH gamma: %s" % gamma)
                    print("RUNNING SIMULATION WITH T: %s" % T)
                    startTime = datetime.datetime.utcnow()
                    results = sim.execute(alpha, beta, T, gamma, time, timestep)
                    print("SIMULATION DONE")
                    endTime = datetime.datetime.utcnow()
                    print("INSERTING INTO DATABASE")
                    manager.insert([
                        {
                            'alpha': alpha,
                            'beta': beta,
                            'T': T,
                            'gamma': gamma,
                            'time': time,
                            'timestep': timestep,
                            'startTime': startTime,
                            'endTime': endTime,
                            'measuredVariable' : dataKey,
                            'data': dataValues
                        } for dataKey, dataValues in results.items()]
                    )

# after all sims have been run
manager.closeConnection()
endClockTime = Time.clock()
print("Controller took: %ss to run" % (endClockTime - startClockTime))