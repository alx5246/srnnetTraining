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

	def uploadData(self, dataFilePath, dbParams, collectionName, urlComponents):
		# upload file to file server
		self._filerserverConnection.upload(dataFilePath, urlParams=urlComponents)

		# add entry to the database
		self._dbManager.openCollection(collectionName)
		dbParams["fileName"] = dataFilePath.split("\\")[-1].split(".")[0]
		dbParams["fileType"] = ".tar.gz"
		dbParams["fileserver_url"] = HTTPRequest.urljoin(*urlComponents)
		self._dbManager.insert(dbParams, insertOne=True)

		# delete local copy
		os.remove(dataFilePath)

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
	def is_within_directory(directory, target):
		
		abs_directory = os.path.abspath(directory)
		abs_target = os.path.abspath(target)
	
		prefix = os.path.commonprefix([abs_directory, abs_target])
		
		return prefix == abs_directory
	
	def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
	
		for member in tar.getmembers():
			member_path = os.path.join(path, member.name)
			if not is_within_directory(path, member_path):
				raise Exception("Attempted Path Traversal in Tar File")
	
		tar.extractall(path, members, numeric_owner=numeric_owner) 
		
	
	safe_extract(tarFile, self._currentDir)

			# copy file to dict
			dirDict = self.makeDirectoryToDictRecursively(os.path.join(self._currentDir, fileName))

			dataDict = extractionFunction(dirDict)
			shutil.rmtree(os.path.join(self._currentDir, fileName))
			outputData.append(dataDict)

		return outputData
