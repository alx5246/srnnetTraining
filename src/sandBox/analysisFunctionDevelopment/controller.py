import datetime
import os
import sys
import multiprocessing


currentDir = os.path.dirname(os.path.realpath(__file__)) # returns a string
layerStripping = [".." for x in range(0, 4)] # a list of ".." strings ["..", "..", "..", ".."]

# RELATIVE TO THE DIRECTORY controller.py is in...currentDir is the string form
# of the file path to the directory that controller.py is in...modify it
dbPath = os.path.join(currentDir, *layerStripping) # equivalent of os.path.join(currentDir, "..", "..", "..", "..")

# let python know where to search
sys.path.extend([currentDir, dbPath])

import dbfileserver as db
import simulation
import time as Time


def wrapList(conn): #, lock):
    fileBaseName = "carlsonOneNeuron_"
    manager = db.DataManager(databaseName="carlsonOneNeuron")
    paramList, proccessNumber = conn.recv()
    conn.close()
    # lock.acquire()
    global currentDir
    brianDir = os.path.join(currentDir, "process_%s_brian_dir" % proccessNumber)
    sim = simulation.SimObject(brianInternalDir=brianDir)
    # lock.release()

    for paramPack, fileCounter in paramList:
        dbParamsDict = dict(paramPack)
        paramPack["iteration"] = fileCounter
        startTime = datetime.datetime.utcnow()
        formattedTime = str(startTime).replace("-", "_").replace(":", "_").replace(" ", "_")
        formattedName =  fileBaseName + str(paramPack["iteration"]) + "_" +\
                         formattedTime
        paramPack["fileName"] = formattedName

        # time and run the simulation
        resultTarGZPath = sim.execute(**paramPack)
        endTime = datetime.datetime.utcnow()

        #
        dbParamsDict['startTime'] = str(startTime)
        dbParamsDict['endTime'] = str(endTime)
        manager.uploadData(resultTarGZPath, dbParamsDict, "STDParameterRange",
                           ["carlsonOneNeuron", "STDParameterRange"])

if __name__ == "__main__":
    startClockTime = Time.clock()

    manager = db.DataManager(databaseName="carlsonOneNeuron")
# manager.openCollection("STDParameterRange")

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



    paramList = []
    allCombinationsFlatList = []

    param1 = {
        "alpha": [0,.1,.2,.3,.4,.5,.6,.7,.8,.9,10],
    }

    param2 = {
        "beta": [0,1,2,4,6,8,10],
    }
    param3 = {
        "gamma": [0,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100],
    }
    param4 = {
        "T": [0,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100],
    }

# print("MONGODB_URI: %s" % os.environ.get("MONGODB_URI"))

    startTime = None
    endTime = None
    results = None

# Make sure to construct the network!! Only make ONE of these objects
    # sim = simulation.SimObject()

    fileBaseName = "carlsonOneNeuron_"
    fileCounter = 1
    for alpha in param1["alpha"]:
        for beta in param2["beta"]:
            for gamma in param3["gamma"]:
                for T in param4["T"]:
                    for simRep in range(10):
                        allCombinationsFlatList.append(
                            ({
                                "alpha": alpha,
                                "beta": beta,
                                "gamma": gamma,
                                "T": T,
                                "time": time,
                                "timestep": timestep,
                            }, fileCounter)
                        )
                        fileCounter += 1

    for i in range(0, len(allCombinationsFlatList), 100000):
        if len(allCombinationsFlatList) - i < 100000:
            paramList.append(allCombinationsFlatList[i:])
        else:
            paramList.append(allCombinationsFlatList[i : i + 100000])

    # print out data to check partitioning
    print("len(paramList): %s" % len(paramList))
    for i in range(0, len(paramList)):
        print("\tlen(paramList[%s]: %s" % (i, len(paramList[i])))

    processes = []
    # lock = multiprocessing.Lock()
    for index in range(0, len(paramList)):
        parent_conn, child_conn = multiprocessing.Pipe()
        processes.append(multiprocessing.Process(target=wrapList,
                                                 # args=(child_conn,lock)))
                                                 args=(child_conn,)))
    # p.map(wrapList, paramList[0][0])
        processes[index].start()
        parent_conn.send((paramList[index], index))

    for index in range(0, len(processes)):
        processes[index].join()

# after all sims have been run
# manager.closeConnection()
    endClockTime = Time.clock()
    print("Controller took: %ss to run" % (endClockTime - startClockTime))