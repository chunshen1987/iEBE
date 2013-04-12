#!/usr/bin/env python
"""
    Use this script to pass a query into a SQLite database.
"""


from sys import argv, stdout, exit

try:
    databaseFilename = argv[1]
    SQLQueryCommandString = " ".join(argv[2:])
    if not SQLQueryCommandString.strip(): raise ValueError("The SQL query should not be empty.")
except:
    print('Usage: databaseQuery.py database_filename "SQL_query"')
    exit()

import sqlite3 as sq3

queryResultCursor = sq3.connect(databaseFilename).execute(SQLQueryCommandString)
stdout.writelines("\t".join(str(element) for element in item) + "\n" for item in queryResultCursor)
