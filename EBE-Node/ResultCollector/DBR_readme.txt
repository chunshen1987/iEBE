
====================================
    Document for the DBR module
====================================

The DBR module is a layer on top of the sqlite3 module, which allows easy to use but less flexible read/write operations to databases.

The main class in the DBR module is the SqliteDB class, which can be used to connect to a database and perform operations. All operations are eventually converted to SQL command, then executed via the sqlite3.execute funtion.

Most of the member function try to hide the details of how to write a correct SQL query from the user; for some functions (e.g. createTableIfNotExists) that does allow the user to directly supply part of a SQL command, an OperationalError exception will be raised if the SQL query is illegal. Another type of exception is the SqliteDBError exception, which is a subclass of the sqlite3.OperationalError, and it is used to represent errors that are not caused by an illegal SQL command from SqliteDB class.

-----------------------------
1. Registering a database
-----------------------------

The module is designed not to constantly keep a connection to preserve resources and a connection is only opened before operations and will be closed immediately afterwards. Since each open operation to the database may potentially fail, the SqliteDB class will not check the validity of the database upon initilization, but only "register" the name of the database. The registered name will be used to any following database operations. To create a SqliteDB object and register a database filename at the same time just pass the filename to the constructor, and the registered database name can be echoed vis the getRegisteredDatabase member function, for example:

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
>>> db.createTableIfNotExists("employee", ("id", "name"), ("int primary key", "char(20)")) # doctest: +ELLIPSIS
<sqlite3.Cursor object ...>
>>> db.doesTableExist("employee")
True
>>> db.doesTableExist("employee1")
False
>>> db.getAllTableNames() == ['employee']
True

The createTableIfNotExists function accepts a table name, a list of column names, and a list of column data types. Inside the type string any valid SQL commands can be used (otherwise an OperationalError exception will be raised). For example, when creating tables, the table name cannot contain blank spaces:
>>> db.createTableIfNotExists("id and name", ("id", "name"), ("int primary key", "char(20)"))
Traceback (most recent call last):
    ...
OperationalError: SQL table name cannot contain blanks

For another example an exception will be raised when there a SQL syntax error:
>>> db.createTableIfNotExists("tableWillNotBeCreated", ("id"), ("int primary keys")) # should be "key", not "keys"
Traceback (most recent call last):
    ...
OperationalError: near "keys": syntax error

The createTableIfNotExists function, as its name suggests, will only create a table when it does not exist:
>>> db.createTableIfNotExists("employee", ("id", "name"), ("int primary key", "char(20)"))
False

----------------------------
3. Read and write tables
----------------------------

Writing to tables is done via the the member function insertIntoTable, named after the name of the corresponding SQL command. It will write a list of values which is its 2nd argument, to the table with the name given by its 1st argument. For example:
>>> db.insertIntoTable("employee", [(1, "Alice"), (2, "Bob"), (3, "Cauchy")]) # doctest: +ELLIPSIS
<sqlite3.Cursor object ...>

It will generate an OperationError if given value does not match the type of the table:
>>> db.insertIntoTable("employee", [(4,), (5,), (6,)])
Traceback (most recent call last):
    ...
OperationalError: table employee has 2 columns but 1 values were supplied

If the given table has not been created yet, the insertIntoTable function provide the possibility to first create the table then perform data insertion, but which requires two additional arguments: a list of column names, and a list of data types; they are exactly the same as those used in createTableIfNotExists function. For example:
>>> db.insertIntoTable("numerical", [(1.3,), (2.4,), (3.7,)], ["value"], ["float"]) # doctest: +ELLIPSIS
<sqlite3.Cursor object ...>

However if the table does not exist and name/type list is not given, it will raise a SqliteDBError exception:
>>> db.insertIntoTable("IDONOTEXIST", [(1.3,), (2.4,), (3.7,)])
Traceback (most recent call last):
    ...
SqliteDBError: table does not exist and cannot be created because name/type list not given

Reading from table is via the member function selectFromTable. Its 1st argument is the name of the table, and its 2nd argument is a list of name of columns (fields) to be read from the table, which defaults to "*" (all fields). It also has an optional named whereClause argument allows to pass any valid SQL where clause (without the keyword "where") to the query. The following are a few examples:

>>> db.selectFromTable("employee")
[(1, u'Alice'), (2, u'Bob'), (3, u'Cauchy')]

>>> db.selectFromTable("employee", "name")
[(u'Alice',), (u'Bob',), (u'Cauchy',)]

>>> db.selectFromTable("employee", whereClause="id>=2")
[(2, u'Bob'), (3, u'Cauchy')]

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
[(u'Alice', 2), (u'Bob', 3), (u'Cauchy', 4)]

For a more fancy example, the following code return name and id for all the employees, but in such a way that their id's are shifted up by the max id in the existing database:
>>> max_id = db.selectFromTable("employee", "max(id)")[0][0]
>>> db.selectFromTable("employee", ["name", "id+%d" % max_id])
[(u'Alice', 4), (u'Bob', 5), (u'Cauchy', 6)]

To obtain a list of all the fields in a table, use the getTableInfo function. It returns a list of tuples of the form (field name, field type). For example:
>>> db.getTableInfo("employee")
[(u'id', u'int'), (u'name', u'char(20)')]

For convenience another function unpackDatabase is also provided. It role is to write out all the tables from a database into separated files. Each file assumes the name of the table it contains; other details for tunable via arguments, like the string used to separate data, whether to include a header, etc. For example:
>>> db.unpackDatabase(sep=",", ext=".dat")
>>> open("employee.dat").readlines()
['#id,name\n', '1,Alice\n', '2,Bob\n', '3,Cauchy\n']

---------------
4. Deletion
---------------

Deletion should always be proceeded by caution. In a real sqlite system there is a command "commit" to actually write changes to the physical database. For the SqliteDB class all operations are automatically followed by a "commit" statement therefore there is no turning back once a disastrous delete operation is issued.

Because of the danger, the SqliteDB class has an internal lock _deletionSafety that is initialized to True. This lock has can be either released or overriden, otherwise a DeletionLockedError will be raised. This lock will be enforced everytime a new database is registered.

The deletion of a table is done through the dropTable function, the following is an example for deletion without the lock:
>>> db._deletionSafety = False # release the safety lock
>>> db.dropTable("employee") # perform deletion
True
>>> db.getAllTableNames() == ['numerical'] # table "employee" is gone
True
>>> db.dropTable("employee") # delete a non-existing table will return false
False

Next is an example of deletion with the lock:
>>> db._deletionSafety = True # relock it
>>> db.dropTable("numerical") # cannot delete a table anymore
Traceback (most recent call last):
    ...
SqliteDBError: deletion prevented by the internal lock

>>> db.dropTable("numerical", overrideSafety=True) # override the safety lock
True

Finally, to delete a database simply delete the file that is "registered" with it, or by calling the deleteDatabase function:
>>> db.deleteDatabase(overrideSafety=True)
True



The END
