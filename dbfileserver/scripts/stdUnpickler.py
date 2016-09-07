# SYSTEM IMPORTS
import pickle

# PYTHON PROJECT IMPORTS


def unpickleKeys(argsList, keysToUnpickle):
    returnList = []
    for elementDict in argsList:
        unpickledDict = dict()
        for key in keysToUnpickle:
            if key in elementDict:
                unpickledDict[key] = pickle.loads(str(elementDict[key]))
        returnList.append(unpickledDict)
    return returnList

def unpickleAll(argsList):
    returnList = []
    for elementDict in argsList:
        unpickledDict = {}
        for key, value in elementDict.items():
            unpickledDict[key] = pickle.loads(str(value))
        returnList.append(unpickledDict)
    return returnList

