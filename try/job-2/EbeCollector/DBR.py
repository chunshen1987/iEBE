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
        self._dbCon = None # reference to the database connection
        self.registerDatabase(fileName)

    def getRegisteredDatabase(self):
        """
            Return the registered database filename.
        """
        return self._registeredDatabase

    def registerDatabase(self, fileName):
        """
            Register and open a database file with "filename".
        """
        self.closeConnection()
        self._registeredDatabase = fileName # used to store the name of the registered database
        self._openConnection()

    def _openConnection(self):
        """
            Open a connection to a database file.
        """
        # check if database is registerred
        if not self._registeredDatabase:
            raise self.SqliteDBError("database not registerd")
        # check if database already open
        if not self._dbCon:
            self._dbCon = sqlite3.connect(self._registeredDatabase)

    def closeConnection(self, discardChanges=False):
        """
            Close the current connection. All modification will be saved upon
            closing, unless "discardChanges" is set to True.
        """
        if self._dbCon:
            if not discardChanges:
                self._dbCon.commit() # write to disk only upon closure
            self._dbCon.close()
            self._dbCon = None

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
        return [item[1] for item in self._executeSQL("select * from sqlite_master") if item[0]=="table"] # returned item: (type, name, tbl_name, rootpage, sql)

    def getTableInfo(self, tableName):
        """
            Return a list of the form ('field', 'type') for all fields in the
            table "tableName".
        """
        return [ (item[1],item[2]) for item in self._executeSQL("pragma table_info(%s)" % tableName) ]

    def doesTableExist(self, tableName):
        """
            Returns True if the table with name "tableName" exists in the
            database.
        """
        return tableName in self.getAllTableNames()

    def createTableIfNotExists(self, tableName, nameAndTypeList):
        """
            Create a table with name "tableName" if it does not exist already.
            The argument "nameAndTypeList" is a list of pair of strings that
            specifies the names and data type of the columns. For example:
            (("id", "integer"), ("name", "text")).
            Returns False if the table already exists.
        """
        # check if table name is legal
        tableName = tableName.strip()
        if " " in tableName:
            raise sqlite3.OperationalError("SQL table name cannot contain blanks")
        # refine input arguments
        if not ListRNew.isIterable(nameAndTypeList[0]):
            nameAndTypeList = [nameAndTypeList]
        # create the table
        if not self.doesTableExist(tableName):
            return self._executeSQL( "create table %s (%s)" % (tableName, ",".join(map(" ".join, nameAndTypeList))) )
        else:
            return False

    def insertIntoTable(self, tableName, valueList):
        """
            Insert values from "valueList" into the table with name "tableName".
            The inserted values have name list "nameList" and types specified in
            "dataTypeStringList". For example:

            tableName = "test"
            valueList = [(1,), (2,), (3,)]
            nameList = ["id"]
            dataTypeStringList = ["int"]

            The table has to be already created (not checked).
        """
        # make valueList doubly nested
        if not ListRNew.isIterable(valueList):
            valueList = [valueList]
        if not ListRNew.isIterable(valueList[0]):
            valueList = [valueList] # nest the level 1 list to level 2 assuming we are given a single value list

        # perform SQL
        dataLength = len(valueList[0]) # get number of elements in a row
        returnValue = self._executeSQL("insert into %s values (%s)" % (tableName, ",".join("?"*dataLength)), valueList, many=True) # perform insertion
        return returnValue


    def selectFromTable(self, tableName, columnNameList="*", whereClause="", groupByClause="", orderByClause=""):
        """
            Return the specified columns with names given in columnNameList from
            the table with name tableName. The columnNameList will be joined
            with a space to be inserted into the SQL query command. The
            whereClause string argument is appended to the query after the
            keyword "where"; the orderByClause string argument is appended to
            the query after the keyword "order by".
        """
        if not ListRNew.isIterable(columnNameList):
            columnNameList = [columnNameList]
        columnNameListSQLString = ",".join(columnNameList) # make is SQL-like
        # perform SQL
        sqlCommand = "select %s from %s" % (columnNameListSQLString, tableName)
        if whereClause:
            sqlCommand += " where " + whereClause
        if groupByClause:
            sqlCommand += " group by " + groupByClause
        if orderByClause:
            sqlCommand += " order by " + orderByClause
        returnValue = self._executeSQL(sqlCommand).fetchall()
        return returnValue

    def dropTable(self, tableName):
        """
            Delete the table with name "tableName". Return True upon success; return False
            if table does not exist.
        """
        # check table existence
        if tableName not in self.getAllTableNames():
            return False
        # perform SQL
        self._executeSQL("drop table %s" % tableName)
        return True

    def deleteDatabase(self, confirmation=False):
        """
            Delete the file associated to the currently registered database. The
            argument confirmation needs to be set to True in order for deletion
            to actually happen.

            Return False if the database file is not found or database is not
            registered; return True if the deletion is completed normally.
        """
        # check safety lock
        if not confirmation:
            raise self.SqliteDBError("deletion not confirmed")
        # check if the database is registered
        db = self.getRegisteredDatabase()
        if not db: return False
        # try to delete it
        if not path.exists(db): return False
        unlink(db)
        return True

    def unpackDatabase(self, sep="\t", writeToFolder=".", ext=".dat", writeHeader=(True, "# ")):
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
