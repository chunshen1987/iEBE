#!/usr/bin/env python
"""

This package serves as another layer between sqlite3 package and the user layer
so that the user need not to know how the sqlite3 works when using database.

For more details on sqlite see:
http://www.sqlite.org/sqlite.html
http://docs.python.org/2/library/sqlite3.html

"""

import sqlite3

import ListRNew

class SqliteDB(object):
    """
        This class wraps up a sqlite connection and serves as an interface for
        I/O.
    """

    def __init__(self, fileName=":memory:"):
        """
            Register the file with "fileName" as the database file.
        """
        self._registeredDatabase = None
        self._dbCon = False
        self.registerDatabase(fileName)

    def getRegisteredDatabase(self):
        """
            Return the registered database filename.
        """
        return self._registeredDatabase

    def registerDatabase(self, fileName):
        """
            Register a database filename.
        """
        self._registeredDatabase = fileName # used to store the name of the database

    def _openConnection(self):
        """
            Open a connection to a database file. Returns True upon success.
        """
        # check if database is registerred
        if not self._registeredDatabase:
            return False
        # close opened connection to prevent port leakage
        self._closeConnection()
        # try to open a new connection
        try:
            self._dbCon = sqlite3.connect(self._registeredDatabase)
        except sqlite3.OperationalError:
            return False
        return True

        return True # open successful

    def _closeConnection(self):
        """
            Close the current connection. Returns False if no connection was
            opened.
        """
        if self._dbCon:
            self._dbCon.commit()
            self._dbCon.close()
            self._dbCon = None
            return True
        else:
            return False

    def _executeSQL(self, cmdString, parameterTuple=(), many=False):
        """
            Execute an SQLite command "cmdString" with parameters
            "parameterTuple". This function does not open connection since
            frequent open/close connection action with each query is very
            inefficient and SQL commands are best to be run in batch.
            When "many" is set to True the executemany function is called.
        """
        # check if the connection is open
        if not self._dbCon:
            return False
        # execute,and return False if error occurs
        try:
            if not many:
                return self._dbCon.execute(cmdString, parameterTuple)
            else:
                return self._dbCon.executemany(cmdString, parameterTuple)
        except sqlite3.OperationalError:
            return False


    def getAllTableNames(self):
        """
            Return a list of table names from the registered database.
        """
        if self._openConnection():
            return [item[1] for item in self._executeSQL("select * from sqlite_master") if item[0]=="table"] # returned item: (type, name, tbl_name, rootpage, sql)
        else:
            return []

    def tableExist(self, tableName):
        """
            Returns True if the table with name "tableName" exists in the
            database.
        """
        return tableName in self.getAllTableNames()

    def _createTableIfNotExists(self, tableName, nameList, dataTypeStringList):
        """
            Create a table with name "tableName" if it does not exist already.
            The argument "nameList" is a list of strings that specifies the
            names of the columns, and the argument "dataTypeStringList" is a
            list of strings that specifies the corresponding data types. For
            example:
            nameList=["id", "name"], dataTypeStringList=["int", "char(8)"].
            Returns True if the table already exists, otherwise False.
        """
        if not self.tableExist(tableName):
            return self._executeSQL( "create table %s (%s)" % (tableName, ",".join(map(" ".join, zip(nameList, dataTypeStringList)))) )

    def insertIntoTable(self, tableName, valueList, nameList=[], dataTypeStringList=[]):
        """
            Insert values from "valueList" into the table with name "tableName".
            The inserted values have name list "nameList" and types specified in
            "dataTypeStringList". For example:

            tableName = "test"
            valueList = [(1), (2), (3)]
            nameList = ["id"]
            dataTypeStringList = ["int"]

            The table will be created if it does not exist; otherwise the
            existing table with the same name will be used, in which case the
            arguments "nameList" and "dataTypeStringList" can be ignored.

            Return False if insertion failed.
        """
        if not self.tableExist(tableName):
            self._createTableIfNotExists(tableName, nameList, dataTypeStringList)
        if self.tableExist(tableName):

            # make valueList doubly nested
            if not ListRNew.isIterable(valueList):
                valueList = [valueList]
            if not ListRNew.isIterable(valueList[0]):
                valueList = [valueList]

            # perform SQL
            dataLength = len(valueList[0]) # get number of elements in a row
            self._openConnection() # open connection
            returnValue = self._executeSQL("insert into %s values (%s)" % (tableName, ",".join("?"*dataLength)), valueList, many=True) # perform insertion
            self._closeConnection() # close connection
            if returnValue:
                return True
            else:
                return False
        else:
            return False # somehow the table cannot be created or was deleted just before insertion


    def selectFromTable(self, tableName, columnNameList="*"):
        """
            Return the specified columns with names given in columnNameList from
            the table with name tableName. The columnNameList will be joined
            with a space to be inserted into the SQL inquiry command.
        """
        if not ListRNew.isIterable(columnNameList):
            columnNameList = [columnNameList]
        columnNameListSQLString = ",".join(columnNameList) # make is SQL-like
        # perform SQL
        self._openConnection()
        returnValue = self._executeSQL("select %s from %s" % (columnNameListSQLString, tableName)).fetchall()
        self._closeConnection()
        return returnValue


