#!/usr/bin/env python

import EbeCollector
import DBR
db = DBR.SqliteDB("tmp.db")
collector = EbeCollector.EbeCollector()
collector.collectEccentricitiesAndRIntegrals("testData", 1, db)
print(db.selectFromTable("ecc_id_lookup"))
collector.createDatabaseFromEventFolders("testData_newStyle")

from shutil import copy
copy("testData_newStyle/CollectedResults.db", "testData_newStyle/CollectedResults_copy.db")
collector.mergeDatabases(DBR.SqliteDB("testData_newStyle/CollectedResults.db"), DBR.SqliteDB("testData_newStyle/CollectedResults_copy.db"))
