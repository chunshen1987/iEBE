#!/usr/bin/env python

import EbeCollector
import DBR
collector = EbeCollector.EbeCollector()
collector.createDatabaseFromEventFolders("testData_oldStyle", collectMode="fromPureHydro")


