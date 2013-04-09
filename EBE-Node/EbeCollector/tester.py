#!/usr/bin/env python

import EbeCollector
import DBR
collector = EbeCollector.EbeCollector()
#collector.createDatabaseFromEventFolders("testData_newStyle", collectMode="fromUrQMD")

db = DBR.SqliteDB("tmp.db")
collector.collectFLowsAndMultiplicities_urqmdBinUtilityFormat("testData", 1, db, multiplicityFactor=0.1)
db.selectFromTable("inte_vn", ("vn_real", "vn_imag"), whereClause="n=3")
