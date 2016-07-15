
# PYTHON PROJECT IMPORTS
import pymongo
import os

# PYTHON PROJECT IMPORTS


class DBManager(object):
    def __init__(self, dbUri=None, databaseName=None, collectionName=None):
        if dbUri is None:
            self.dbUri = os.environ.get("MONGODB_URI")
            if self.dbUri is None:
                raise Exception("Could not connect to mongo db. MONGODB_URI env var not set")
        else:
            self.dbUri = dbUri
        self.client = pymongo.MongoClient(self.dbUri)
        if databaseName is None:
            self.databaseName = os.environ.get("MONGODB_NAME")
            if self.databaseName is None:
                raise Exception("Cannot get database from server. MONGODB_NAME env var not set")
        else:
            self.databaseName = databaseName
        self.database = self.client[self.databaseName]
        self.collection = self.database[collectionName] if collectionName is not None else None
        self.collectionName = collectionName

    def openCollection(self, collectionName):
        if collectionName != self.collectionName:
            self.collection = self.database[collectionName]
            self.collectionName = collectionName

    def query(self, queryParams, returnOne=False, sortScheme=None, hook=None):
        results = []
        if returnOne:
            results = [self.collection.find_one(queryParams)]
        else:
            if sortScheme is not None:
                results = [x for x in self.collection.find(queryParams).sort(sortScheme)]
            else:
                results = [x for x in self.collection.find(queryParams)]
        if hook is not None:
            return hook(results)
        return results

    def insert(self, data, insertOne=False):
        if insertOne:
            self.collection.insert_one(data)
        else:
            self.collection.insert_many(data)

    def update(self, queryParams, updateData, updateOne=False):
        if updateOne:
            self.collection.update_one(queryParams, updateData)
        else:
            self.collection.update_many(queryParams, updateData)

    def getAllKeysInCollection(self, paramsToIgnore=[]):
        allParams = set()
        for doc in self.query({}):
            for key in doc.keys():
                if key not in allParams and key not in paramsToIgnore:
                    allParams.add(key)
        return allParams

    def getAllKeysAndValuesInCollection(self, paramsToIgnore=[]):
        allParamsDict = {}
        for key in self.getAllKeysInCollection(paramsToIgnore):
            allParamsDict[key] = self.collection.distinct(key)
        return allParamsDict

    def delete(self, queryParams, deleteOne=False):
        if deleteOne:
            self.collection.delete_one(queryParams)
        else:
            self.collection.delete_many(queryParams)

    def closeConnection(self):
        self.database = None
        self.collection = None
        self.client.close()
