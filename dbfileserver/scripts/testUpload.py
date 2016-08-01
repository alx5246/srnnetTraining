import os
import sys

currentDir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(currentDir)
import DataManager


x = DataManager.DataManager(databaseName="testDB")
fileserverPath = os.path.realpath(os.path.join(currentDir, ".."))
dataDict = {
	"init": os.path.join(fileserverPath, "__init__.py"),
	"scripts": {
		"HTTPRequest": os.path.join(currentDir, "HTTPRequest.py"),
		"DBManager": os.path.join(currentDir, "DBManager.py")
	}
}
x.uploadData(dataDict, {"upload": 1}, "testCollection", "testUpload", ["carlsonOneNeuron", "STDParameterRange"])

def extractionFunction(dirDict):
	return dirDict

def prettyPrintDict(dictionary, offset):
	print(offset + "{")
	offset = "\t"
	for key in dictionary.keys():
		if isinstance(dictionary[key], dict):
			print(offset + key + ":")
			prettyPrintDict(dictionary[key], offset + "\t")
		else:
			print(offset + key + ": " + dictionary[key])

dataList = x.downloadData("testCollection", {"upload": 1}, extractionFunction)
print("len(dataList): %s" % len(dataList))
# for data in dataList:
# 	prettyPrintDict(data, "")
