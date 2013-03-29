#!/usr/bin/env python

import EbeCollector
import DBR
db = DBR.SqliteDB("tmp.db")
collector = EbeCollector.EbeCollector()
collector.collectEccentricitiesAndRIntegrals("testData", 1, db)
print(db.selectFromTable("ecc_id_lookup"))
collector.collectFLowsAndMultiplicities_binUtilityFormat("testData", 1, db, multiplicityFactor=0.1)
