from EbeCollector import EbeDBReader


# The following strings are used as query parameters.
EXPR = "expr"
DATABASE = "database"


class QueryBridge(object):
    def __init__(self, database_name=""):
        self.database_name = None
        self.reader = None

        self.set_datebase(database_name)

    def set_database(self, database_name):
        if database_name != self.database_name:
            self.database_name = database_name
            self.reader = EbeDBReader(database_name)

    def evaluate_expression(self, expr):
        return self.reader.evaluateExpressionOnly(expr)
