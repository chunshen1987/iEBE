#!/usr/bin/env python
"""

This is the test file goes with the DBR package. These test functions are meaned
to be manually called or via the PyTestFunctionCaller module. Note that the
execution of later tests may depends on the success of the execution of former
function.

"""


from os import path, unlink
import DBR

def test_01_SqliteDB_constractor():
    test_file = "test_tmp.db"
    if path.exists(test_file): unlink(test_file)
    testDB = DBR.SqliteDB(test_file)
    assert testDB.getRegisteredDatabase() == test_file
    return True

def test_02_SqliteDB__openAndCloseConnection():
    test_file = "test_tmp.db" # valid filename
    testDB = DBR.SqliteDB(test_file)
    assert testDB._openConnection() == True
    assert testDB._closeConnection() == True
    unlink(test_file)
    test_file = path.join(*(["tmp"]*8)) # invalid filename
    testDB = DBR.SqliteDB(test_file)
    assert testDB._openConnection() == False
    assert testDB._closeConnection() == False
    return True

def test_03_SqliteDB__executeSQL():
    test_file = "test_tmp.db"
    testDB = DBR.SqliteDB(test_file)
    assert testDB._executeSQL("abcde") == False # invalid command
    assert testDB._executeSQL("create table test (id int)") == False # valid command but connection not open
    testDB._openConnection()
    assert testDB._executeSQL("create table test (id int)") != False # valid command with open connection
    assert testDB._executeSQL("insert into test values (3)") != False
    assert testDB._executeSQL("insert into test values (?)", [(4,),(5,),(6,),(7,)], many=True) != False
    testDB._closeConnection()
    unlink(test_file)
    return True

def test_04_SqliteDB_getAllTableNames():
    test_file = "test_tmp.db"
    testDB = DBR.SqliteDB(test_file)
    assert testDB.getAllTableNames() == [] # none since the database is empty
    testDB._openConnection()
    testDB._executeSQL("create table test1 (id int)")
    testDB._executeSQL("create table test2 (id int)")
    assert testDB.getAllTableNames() == ["test1", "test2"] # now two databases
    testDB._closeConnection()
    unlink(test_file)
    return True

def test_05_SqliteDB__createTableIfNotExists():
    test_file = "test_tmp.db"
    testDB = DBR.SqliteDB(test_file)
    testDB._createTableIfNotExists("test", ["id"], ["int"]) # create a new table
    assert "test" in testDB.getAllTableNames() # it should contain the test table now
    testDB._createTableIfNotExists("test", ["id"], ["int"]) # creating the same table should have no effects
    assert testDB.getAllTableNames() == ["test"] # should have only 1 table
    testDB._closeConnection()
    unlink(test_file)
    return True

def test_06_SqliteDB_insertIntoTable():
    test_file = "test_tmp.db"
    testDB = DBR.SqliteDB(test_file)
    assert not testDB._dbCon # the database should be closed
    assert testDB.insertIntoTable("test", [1], ["id"], ["int"]) # single value insertion
    assert testDB.insertIntoTable("test", [(1,), (2,)], ["id"], ["int"]) # multiple value insertion
    assert testDB.insertIntoTable("test", [1, 2], ["id"], ["int"]) == False # wrong syntax, should fail
    assert testDB.insertIntoTable("test1", [[1,"a"], [2, "b"]], ["id", "name"], ["int", "char(8)"]) # this complicated command should work
    assert not testDB._dbCon # the database should be closed
    unlink(test_file)
    return True

def test_07_SqliteDB_readFromTable():
    test_file = "test_tmp.db"
    testDB = DBR.SqliteDB(test_file)
    testDB.insertIntoTable("test", [[1,"a"], [2, "b"]], ["id", "name"], ["int", "char(8)"]) # this complicated command should work
    assert testDB.selectFromTable("test", ["id", "name"]) == [(1,"a"), (2, "b")] # what you save is what you read
    assert testDB.selectFromTable("test") == [(1,"a"), (2, "b")] # test the default value
    assert testDB.selectFromTable("test", "id") == [(1,), (2,)] # check selection from single column
    unlink(test_file)
    return True

if __name__ == '__main__':
    import PyTestFunctionCaller
    PyTestFunctionCaller.callTestFunctionsInModule("DBR_test", interceptException=False)
