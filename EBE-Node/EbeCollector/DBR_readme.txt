
====================================
    Document for the DBR module
====================================

The DBR module is a layer on top of the sqlite3 module, which allows easy to use but less flexible read/write operations to databases.

The main class in the DBR module is the SqliteDB class, which can be used to connect to a database and perform operations. All operations are eventually converted to SQL command, then executed via the sqlite3.execute funtion.

Most of the member function try to hide the details of how to write a correct SQL query from the user; for some functions (e.g. createTableIfNotExists) that does allow the user to directly supply part of a SQL command, an OperationalError exception will be raised if the SQL query is illegal. Another type of exception is the SqliteDBError exception, which is a subclass of the sqlite3.OperationalError, and it is used to represent errors that are not caused by an illegal SQL command from SqliteDB class.

-----------------------------
1. Registering a database
-----------------------------

The module is designed to constantly keep a connection open for efficient data reading/writing. A connection is open when a database file is "registered", which can be done during initialization or by calling the registerDatabase function. Once a database file is registered, it will keep open until it is been closed explicitly, or when another registration is issued (by the registerDatabase function). If a connection is closed, the database file is still registered and any subsequenctial read/write operations will re-open the connection. Changes are only written to the database file upon closing so a close action *MUST* be performed (either via closeConnection function or registerDatabase function) in order to make changes to the real database files. See the section "commit changes" for more information.

To create a SqliteDB object and register a database filename at the same time just pass the filename to the constructor, and the registered database name can be echoed vis the getRegisteredDatabase member function, for example:
>>> from DBR import *
>>> db = SqliteDB("myDatabase.db")
>>> db.getRegisteredDatabase() == "myDatabase.db"
True

When absent, the argument in the constructor will defaults to ":memory:" which is a temporary database in memory that is idea for performing data analysis without modifying the physical database.

The registered database name can be changed at any time via the member function registerDatabase, for example:

>>> db.registerDatabase("myAnotherDatabase.db")
>>> db.getRegisteredDatabase() == "myAnotherDatabase.db"
True

---------------------------
2. Operations on tables
---------------------------

To create a table with a given name use the createTableIfNotExists function, to check if a table exists use the doesTableExist function, and to return a list of all the table names use the getAllTableNames member function. For example:

>>> db.registerDatabase("tmp_test_database.db")
>>> db.createTableIfNotExists("employee", (("id", "integer"), ("name", "text"))) # doctest: +ELLIPSIS
<sqlite3.Cursor object ...>
>>> db.doesTableExist("employee")
True
>>> db.doesTableExist("employee1")
False
>>> db.getAllTableNames() == ['employee']
True

The createTableIfNotExists function accepts a table name, a field name-type list. Both strings must be valid SQL commands, otherwise an OperationalError exception will be raised. For example, when creating tables, the table name cannot contain blank spaces:
>>> db.createTableIfNotExists("id and name", (("id", "integer"), ("name", "text")))
Traceback (most recent call last):
    ...
OperationalError: SQL table name cannot contain blanks

For another example an exception will be raised when there a SQL syntax error:
>>> db.createTableIfNotExists("tableWillNotBeCreated", ("id", "int primary keys")) # should be "key", not "keys"
Traceback (most recent call last):
    ...
OperationalError: near "keys": syntax error

The createTableIfNotExists function, as its name suggests, will only create a table when it does not exist:
>>> db.createTableIfNotExists("employee", (("id", "integer"), ("name", "text")))
False

----------------------------
3. Read and write tables
----------------------------

Writing to tables is done via the the member function insertIntoTable, named after the name of the corresponding SQL command. It will write a list of values which is its 2nd argument, to the table with the name given by its 1st argument. For example:
>>> db.insertIntoTable("employee", [(1, "Alice"), (2, "Bob"), (3, "Cauchy")]) # doctest: +ELLIPSIS
<sqlite3.Cursor object ...>

When writing only one value entry (row) into the database, the value list does not need to be doubly nested, for example:
>>> db.insertIntoTable("employee", (4, "David")) # doctest: +ELLIPSIS
<sqlite3.Cursor object ...>

It will generate an OperationError if given value does not match the type of the table:
>>> db.insertIntoTable("employee", [(4,), (5,), (6,)])
Traceback (most recent call last):
    ...
OperationalError: table employee has 2 columns but 1 values were supplied

If the given table has not been created yet, the insertIntoTable function raise an exception. For example:
>>> db.insertIntoTable("numerical", [(1.3,), (2.4,), (3.7,)])
Traceback (most recent call last):
    ...
OperationalError: no such table: numerical

Reading from table is via the member function selectFromTable. Its 1st argument is the name of the table, and its 2nd argument is a list of name of columns (fields) to be read from the table, which defaults to "*" (all fields). It also has an optional named whereClause argument allows to pass any valid SQL where clause (without the keyword "where") to the query. The following are a few examples:

>>> db.selectFromTable("employee")
[(1, u'Alice'), (2, u'Bob'), (3, u'Cauchy'), (4, u'David')]

>>> db.selectFromTable("employee", "name")
[(u'Alice',), (u'Bob',), (u'Cauchy',), (u'David',)]

>>> db.selectFromTable("employee", whereClause="id>=2")
[(2, u'Bob'), (3, u'Cauchy'), (4, u'David')]

All syntax errors with be feedbacked via exceptions:

>>> db.selectFromTable("employee", "age") # no such field as "age"
Traceback (most recent call last):
    ...
OperationalError: no such column: age

>>> db.selectFromTable("employer") # no such table
Traceback (most recent call last):
    ...
OperationalError: no such table: employer

>>> db.selectFromTable("employer", whereClause="hello world!") # illegal where clause
Traceback (most recent call last):
    ...
OperationalError: near "world": syntax error

Any valid SQL syntax will be passed along to the underlayer, for example, to shift up all the employees' id by 1 we can do:
>>> db.selectFromTable("employee", ["name", "id+1"])
[(u'Alice', 2), (u'Bob', 3), (u'Cauchy', 4), (u'David', 5)]

For a more fancy example, the following code return name and id for all the employees, but in such a way that their id's are shifted up by the max id in the existing database:
>>> max_id = db.selectFromTable("employee", "max(id)")[0][0]
>>> db.selectFromTable("employee", ["name", "id+%d" % max_id])
[(u'Alice', 5), (u'Bob', 6), (u'Cauchy', 7), (u'David', 8)]

To obtain a list of all the fields in a table, use the getTableInfo function. It returns a list of tuples of the form (field name, field type). For example:
>>> db.getTableInfo("employee")
[(u'id', u'integer'), (u'name', u'text')]

For convenience another function unpackDatabase is also provided. It role is to write out all the tables from a database into separated files. Each file assumes the name of the table it contains; other details for tunable via arguments, like the string used to separate data, whether to include a header, etc. For example:
>>> db.unpackDatabase(sep=",", ext=".dat")
>>> open("employee.dat").readlines()
['# id,name\n', '1,Alice\n', '2,Bob\n', '3,Cauchy\n', '4,David\n']

---------------
4. Deletion
---------------

In a real sqlite system there is a command "commit" to actually write changes to the physical database. For the SqliteDB class, all operations are written to the database upon closing (however, see the section "commit changes") and there is no explicit commit operation.

The deletion of a table is done through the dropTable function, the following is an example:
>>> db.dropTable("employee")
True
>>> db.getAllTableNames() == [] # table "employee" is gone
True
>>> db.dropTable("employee") # delete a non-existing table will return false
False

Finally, to delete a database simply delete the file that is "registered" with it, or by calling the deleteDatabase function, which needs a named "confirmation" argument as a warning.
>>> db.deleteDatabase(confirmation=True)
True

---------------------
5. Commit changes
---------------------

As mentioned upon closing of a connection all changes will ge written to the database file, which will happen either when the closeConnection is explicitly called or when another database file is registered. However if changes are not meant to be written to the database the function, the function closeConnection function also accepts a named argument "discardChanges" which, when setting to True, will close the connection without commit changes (however if a table is created, this action is recorded; though any write action follows it will be discarded). The following examples show this behavior.

First is an example that writes to a database before closing.
>>> db.registerDatabase("myDatabase.db") # this opens a connection
>>> db.createTableIfNotExists("integer", ("i", "int")) # creates a single table #doctest: +ELLIPSIS
<sqlite3.Cursor object at ...>
>>> db.insertIntoTable("integer", ((1,), (2,), (3,), (4,), (5,))) # insert values into the table #doctest: +ELLIPSIS
<sqlite3.Cursor object at ...>
>>> db.closeConnection() # changes are written to the database file
>>> db.unpackDatabase() # this action will re-open the connection to the same database file and unpack it
>>> open("integer.dat").readlines()
['# i\n', '1\n', '2\n', '3\n', '4\n', '5\n']

>>> db.deleteDatabase(confirmation=True)
True

Next is an example that discard changes upon closing.
>>> db.registerDatabase("myDatabase.db") # this opens a connection
>>> db.createTableIfNotExists("integer", ("i", "int")) # creates a single table #doctest: +ELLIPSIS
<sqlite3.Cursor object at ...>
>>> db.insertIntoTable("integer", ((1,), (2,), (3,), (4,), (5,))) # insert values into the table #doctest: +ELLIPSIS
<sqlite3.Cursor object at ...>
>>> db.closeConnection(discardChanges=True) # changes discarded
>>> db.unpackDatabase() # this action will re-open the connection to the same database file and unpack it
>>> open("integer.dat").readlines() == ['# i\n'] # only the table is visible, no data
True
>>> db.deleteDatabase(confirmation=True)
True


-------------
Clean ups
-------------
>>> unlink("myAnotherDatabase.db")
>>> unlink("employee.dat")
>>> unlink("integer.dat")

The END
