from os import path
from EbeCollector import EbeDBReader


# The following strings are used as query parameters.
EXPR_PARAM = "expr"
DATABASE_PARAM = "database"

DATABASE_RELATIVE_PATH = "databases"


class QueryBridge(object):
    def __init__(self, database_name=""):
        self.database_name = ""
        self.reader = None

        self.set_datebase(database_name)

    def set_database(self, database_name):
        if not database_name.endwith(".db"):
            database_name += ".db"

        if database_name != self.database_name:
            self.database_name = database_name
            self.reader = EbeDBReader(path.join(DATABASE_RELATIVE_PATH,
                                                database_name))

    def evaluate_expression(self, expr):
        return self.reader.evaluateExpressionOnly(expr)
