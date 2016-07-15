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

manager = db.DBManager(databaseName="carlsonOneNeuron")
manager.openCollection("STDParameterRange")

time = 10   #seocnds
timestep = 100  #In units of microseconds
gamma = 50
beta = 1
T = 5 # units of seconds

#Setup paramter lists
paramsDict = {
    "alpha": [0.4],
}

print("MONGODB_URI: %s" % os.environ.get("MONGODB_URI"))

startTime = None
endTime = None
results = None
for alpha in paramsDict["alpha"]:

    print("RUNNING SIMULATION WITH ALPHA: %s" % alpha)
    startTime = datetime.datetime.utcnow()
    results = simulation.execute(alpha, beta, T, gamma, time, timestep)
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
