#!/usr/bin/env python
"""

This package serves as another layer between sqlite3 package and the user layer
so that the user need not to know how the sqlite3 module works in order to use
database.

Only comments related to developments are included in this file; for document
see DBR_readme.txt.

For more details on sqlite see:
http://www.sqlite.org/sqlite.html
http://docs.python.org/2/library/sqlite3.html

"""

import sqlite3
from os import unlink, path
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
        self._registeredDatabase = None # stores the filename for the database
        self._dbCon = False # reference to the database connection
        self._deletionSafety = True # extra caution will be used during deletion when set to True
        self._counter = 0 # debugging
        self.registerDatabase(fileName)

    def getRegisteredDatabase(self):
        """
            Return the registered database filename.
        """
        return self._registeredDatabase

    def registerDatabase(self, fileName):
        """
            Register a database filename. Note that this action will
            delteion-safety lock the newly registered database, even if the
            newly registered database is the same as the old one.
        """
        self._registeredDatabase = fileName # used to store the name of the database
        self._deletionSafety = True

    def _openConnection(self):
        """
            Open a connection to a database file. Returns True upon success.
        """
        # check if database is registerred
        if not self._registeredDatabase: raise self.SqliteDBError("database not registerd")
        # close opened connection to prevent port leakage
        self._closeConnection()
        # throw exceptions upon errors
        self._dbCon = sqlite3.connect(self._registeredDatabase)
        return True

    def _closeConnection(self):
        """
            Close the current connection. Returns False if no connection was
            opened. All modification will be saved upon closure.
        """
        if self._dbCon:
            self._dbCon.commit()
            self._counter += 1
            print(self._counter)
            self._dbCon.close()
            self._dbCon = None
            return True
        else:
            return False

    def _executeSQL(self, cmdString, parameterTuple=(), many=False):
        """
            Execute an SQLite command "cmdString" with parameters
            "parameterTuple". This function will open a connection if it is not
            already open, but it will not close the connection since frequent
            open/close connection action with each query is very inefficient and
            SQL commands are best to be run in batch. When "many" is set to True
            the executemany function is called.
        """
        # check if the connection is open
        if not self._dbCon:
            self._openConnection()
        # execute,and throw exceptions if error occurs
        try:
            if not many:
                return self._dbCon.execute(cmdString, parameterTuple)
            else:
                return self._dbCon.executemany(cmdString, parameterTuple)
        except sqlite3.OperationalError:
            print("Error executing: %s" % cmdString)
            print("With parameters:")
            print(parameterTuple)
            raise

    def getAllTableNames(self):
        """
            Return a list of table names from the registered database.
        """
        self._openConnection()
        return [item[1] for item in self._executeSQL("select * from sqlite_master") if item[0]=="table"] # returned item: (type, name, tbl_name, rootpage, sql)
        self._closeConnection()

    def getTableInfo(self, tableName):
        """
            Return a list of the form ('field', 'type') for all fields in the
            table "tableName".
        """
        self._openConnection()
        return [ (item[1],item[2]) for item in self._executeSQL("pragma table_info(%s)" % tableName) ]
        self._closeConnection()

    def doesTableExist(self, tableName):
        """
            Returns True if the table with name "tableName" exists in the
            database.
        """
        return tableName in self.getAllTableNames()

    def createTableIfNotExists(self, tableName, nameList, dataTypeStringList):
        """
            Create a table with name "tableName" if it does not exist already.
            The argument "nameList" is a list of strings that specifies the
            names of the columns, and the argument "dataTypeStringList" is a
            list of strings that specifies the corresponding data types. For
            example:
            nameList=["id", "name"], dataTypeStringList=["int", "char(8)"].
            Returns True if the table already exists, otherwise False.
        """
        # check if table name is legal
        tableName = tableName.strip()
        if " " in tableName:
            raise sqlite3.OperationalError("SQL table name cannot contain blanks")
        # refine input arguments
        if not ListRNew.isIterable(nameList):
            nameList = [nameList]
        if not ListRNew.isIterable(dataTypeStringList):
            dataTypeStringList = [dataTypeStringList]
        # create the table
        if not self.doesTableExist(tableName):
            return self._executeSQL( "create table %s (%s)" % (tableName, ",".join(map(" ".join, zip(nameList, dataTypeStringList)))) )
        else:
            return False

    def insertIntoTable(self, tableName, valueList, nameList=None, dataTypeStringList=None):
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
        if not self.doesTableExist(tableName):
            if nameList:
                self.createTableIfNotExists(tableName, nameList, dataTypeStringList)
            else:
                raise self.SqliteDBError("table does not exist and cannot be created because name/type list not given")
        if self.doesTableExist(tableName):

            # make valueList doubly nested
            if not ListRNew.isIterable(valueList):
                valueList = [valueList]
            if not ListRNew.isIterable(valueList[0]):
                valueList = [[x] for x in valueList] # nest the level 1 list to level 2

            # perform SQL
            dataLength = len(valueList[0]) # get number of elements in a row
            self._openConnection() # open connection
            returnValue = self._executeSQL("insert into %s values (%s)" % (tableName, ",".join("?"*dataLength)), valueList, many=True) # perform insertion
            self._closeConnection() # close connection
            return returnValue
        else:
            return False # somehow the table cannot be created or was deleted just before insertion


    def selectFromTable(self, tableName, columnNameList="*", whereClause=""):
        """
            Return the specified columns with names given in columnNameList from
            the table with name tableName. The columnNameList will be joined
            with a space to be inserted into the SQL query command. The
            whereClause string argument is appended to the query after the
            keyword "where".
        """
        if not ListRNew.isIterable(columnNameList):
            columnNameList = [columnNameList]
        columnNameListSQLString = ",".join(columnNameList) # make is SQL-like
        # perform SQL
        self._openConnection()
        sqlCommand = "select %s from %s" % (columnNameListSQLString, tableName)
        if whereClause:
            sqlCommand += " where " + whereClause
        returnValue = self._executeSQL(sqlCommand).fetchall()
        self._closeConnection()
        return returnValue

    def dropTable(self, tableName, overrideSafety=False):
        """
            Delete the table with name "tableName". When self._deletionSafety is
            set to True, the argument "overrideSafty" needs to be set to True in
            order for deletion to happen. Return True upon success; return False
            if table does not exist or deletion is not permitted.
        """
        # check safety lock
        if self._deletionSafety and not overrideSafety:
            raise self.SqliteDBError("deletion prevented by the internal lock")
        # check table existence
        if tableName not in self.getAllTableNames():
            return False
        # perform SQL
        self._openConnection()
        self._executeSQL("drop table %s" % tableName)
        self._closeConnection()
        return True

    def deleteDatabase(self, overrideSafety=False):
        """
            Delete the file associated to the currently registered database.
            When self._deletionSafety is set to True, the argument
            "overrideSafety" needs to be set to True in order for deletion to
            happen.

            Return False if the database file is not found or database is not
            registered; return True if the deletion is completed normally.
        """
        # check safety lock
        if self._deletionSafety and not overrideSafety:
            raise self.SqliteDBError("deletion prevented by the internal lock")
        # check if the database is registered
        db = self.getRegisteredDatabase()
        if not db: return False
        # try to delete it
        if not path.exists(db): return False
        unlink(db)
        return True

    def unpackDatabase(self, sep=" ", writeToFolder=".", ext=".dat", writeHeader=(True, "#")):
        """
            This function writes all content of of a database into files. Each
            table will be written into a single file with table name as the
            filename and "ext" as extension; the data in the files will be
            separated by "sep".

            If writeHeader[0] is set to True then a header containing the names
            of the fiels will be written to the 1st line of the data file, and
            it will be written after a symbol writeHeader[1].
        """
        # get all table names
        tableNames = self.getAllTableNames()

        # write out table one by one
        for aTable in tableNames:
            # create a file with the correpsonding name
            with open(path.join(writeToFolder, aTable+ext), "w") as tableFile:
                # write header
                if writeHeader[0]:
                    fieldNames = [ item[0] for item in self.getTableInfo(aTable) ]
                    tableFile.write(writeHeader[1] + sep.join(fieldNames) + "\n")
                for aRecond in self.selectFromTable(aTable):
                    tableFile.write(sep.join( map(str, aRecond) ) + "\n")

    class SqliteDBError(sqlite3.OperationalError):
        """
            General error exception encountered during database operations.
        """
        pass

if __name__ == '__main__':
    import doctest
    doctest.testfile("DBR_readme.txt")
