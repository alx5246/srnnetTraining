
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
                results = self.collection.find(queryParams).sort(sortScheme)
            else:
                results = self.collection.find(queryParams)
        if hook is not None:
            return hook(results)
        return results

    def insert(self, data, insertOne=False):
        if insertOne:
            self.collection.insert_one(data)
        else:
            self.collection.insert_many(data)

    def closeConnection(self):
        self.database = None
        self.collection = None
        self.client.close()
