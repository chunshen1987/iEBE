from os import path
from EbeCollector import EbeDBReader


# The following strings are used as query parameters.
EXPR_PARAM = "expr"
DATABASE_PARAM = "database"

DATABASE_RELATIVE_PATH = "databases"


class QueryBridge(object):
    def __init__(self):
        self.database_name = ""
        self.reader = None

    def set_database(self, database_name):
        if not database_name.endswith(".db"):
            database_name += ".db"

        if database_name != self.database_name:
            # "str" is used to convert type from unicode to str.
            self.database_name = str(database_name)
            self.reader = EbeDBReader(path.join(DATABASE_RELATIVE_PATH,
                                                self.database_name))

    def evaluate_expression(self, expr):
        if self.reader:
            return self.reader.evaluateExpressionOnly(expr)
        else:
            raise Exception("Invalid database.")
