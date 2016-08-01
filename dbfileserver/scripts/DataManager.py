# SYSTEM IMPORTS
import os
import sys
import shutil
import tarfile


currentDir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(currentDir)


# PYTHON PROJECT IMPORTS
import DBManager
import HTTPRequest


class DataManager(object):
	def __init__(self, dbUri=None, databaseName=None, fileserverUri=None):
		self._filerserverConnection = None
		self._dbManager = DBManager.DBManager(dbUri, databaseName=databaseName)

		if fileserverUri is None:
			self._filerserverConnection = HTTPRequest.HTTPRequest(os.environ["FILESERVER_URI"])
		else:
			self._filerserverConnection = HTTPRequest.HTTPRequest(fileserverUri)
		self._currentDir = os.path.dirname(os.path.realpath(__file__))

	def convertDictToDirRecursively(self, dataDict, fileName, currentDir):
		if not os.path.exists(currentDir):
			os.makedirs(currentDir)
		for dataFile in dataDict.keys():
			# we need to make a subdirectory
			if isinstance(dataDict[dataFile], dict):
				self.convertDictToDirRecursively(dataDict[dataFile], fileName + "_" + dataFile, os.path.join(currentDir, dataFile))
			else:
				# write data to temp file
				# get file type of the file
				filename, file_extension = os.path.splitext(dataDict[dataFile])
				# copy file to directory
				shutil.copyfile(dataDict[dataFile], os.path.join(currentDir, dataFile + file_extension))
				# delete original copy
				# os.remove(dataDict[dataFile])

	def uploadData(self, dataDict, dbParams, collectionName, fileName, urlComponents):
		self.convertDictToDirRecursively(dataDict, fileName, os.path.join(self._currentDir, fileName))

		# create the .tar.gz file
		with tarfile.open(os.path.join(self._currentDir, fileName + ".tar.gz"), "w:gz") as tarFile:
			tarFile.add(os.path.join(self._currentDir, fileName), arcname=fileName)

		# upload file to file server
		self._filerserverConnection.upload(os.path.join(self._currentDir, fileName + ".tar.gz"), urlParams=urlComponents)

		# add entry to the database
		self._dbManager.openCollection(collectionName)
		dbParams["fileName"] = fileName
		dbParams["fileType"] = ".tar.gz"
		dbParams["fileserver_url"] = HTTPRequest.urljoin(*urlComponents)
		self._dbManager.insert(dbParams, insertOne=True)

		# delete local copy
		os.remove(os.path.join(self._currentDir, fileName + ".tar.gz"))
		shutil.rmtree(os.path.join(self._currentDir, fileName))

	def makeDirectoryToDictRecursively(self, dirPath):
		dirDict = {}
		for element in os.listdir(dirPath):
			if os.path.isdir(os.path.join(dirPath, element)):
				dirDict[element] = self.makeDirectoryToDictRecursively(os.path.join(dirPath, element))
			else:
				dirDict[element] = os.path.join(dirPath, element)
		return dirDict

	def downloadData(self, collectionName, queryDict, extractionFunction, sortScheme=None):

		self._dbManager.openCollection(collectionName)
		records = self._dbManager.query(queryDict, sortScheme=sortScheme)
		outputData = []
		print("len(records): %s" % len(records))
		for record in records:
			url = record["fileserver_url"]
			fileName = record["fileName"]
			fileType = record["fileType"]
			urlComponents = url.split("/") + [fileName + fileType]

			# download file
			self._filerserverConnection.download(os.path.join(self._currentDir, fileName + fileType), urlParams=urlComponents)

			# extract file
			with tarfile.open(os.path.join(self._currentDir, fileName + fileType), "r:gz") as tarFile:
				tarFile.extractall(self._currentDir)

			# copy file to dict
			dirDict = self.makeDirectoryToDictRecursively(os.path.join(self._currentDir, fileName))

			dataDict = extractionFunction(dirDict)
			shutil.rmtree(os.path.join(self._currentDir, fileName))
			outputData.append(dataDict)

		return outputData
