#!/usr/bin/env python

from DBR import SqliteDB

fromDB = SqliteDB("/home/qiu/Downloads/Pb_30_40_Glb_More_2000.db")
toDB = SqliteDB("/home/qiu/Downloads/Pb_30_40_Glb_More_2000_truncated.db")

for aTable in fromDB.getAllTableNames():
    print(aTable)
    # first copy table structure
    firstCreation = toDB.createTableIfNotExists(aTable, fromDB.getTableInfo(aTable))
    if "lookup" in aTable: # just copy
        toDB.insertIntoTable(aTable, fromDB.selectFromTable(aTable))
    else: # truncate
        toDB.insertIntoTable(aTable, fromDB.selectFromTable(aTable, whereClause="event_id<=1000"))
toDB.closeConnection() # commit
