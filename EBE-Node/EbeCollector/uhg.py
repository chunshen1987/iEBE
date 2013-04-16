#!/usr/bin/env python
"""
    This module setups a environment for later convenient calculations using the
    EbeDBReader class.
"""

from EbeCollector import EbeDBReader

e = None
c = None

def use(database):
    """
        Create a EbeDBReader object and link the factory functions for
        evaluateExpression and evaluateExpressionOnly as uhg_check and uhg.
    """
    global eval
    global check
    c, e = EbeDBReader(database).getFactoryFunctions()

def help():
    """
        Display a short help message.
    """
    pass


help()

if __name__ == '__main__':
    pass