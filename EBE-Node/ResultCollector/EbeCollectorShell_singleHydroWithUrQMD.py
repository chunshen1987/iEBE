#!/usr/bin/env python
"""
    This is one of the shells to the EbeCollector class. This one
    creates a database using data from subfolders containing multiple
    hybrid (hydro+UrQMD) events.
"""

from sys import argv

if len(argv)<2:
    print("Usage: shell result_folder_name [sub_folder_pattern] [database_name]")
else:
    if len(argv)==2: subfolderPattern = "event-(\d*)"
    from EbeCollector import EbeCollector
    EbeCollector().createDatabaseFromEventFolders()
