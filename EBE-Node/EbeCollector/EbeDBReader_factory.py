#!/usr/bin/env python
"""
    This module setups a environment for later convenient calculations using the
    EbeDBReader class.
"""

from numpy import *
from EbeCollector import EbeDBReader

uhg = None
uhg_check = None

def use(database):
    """
        Create a EbeDBReader object and link the factory functions for
        evaluateExpression and evaluateExpressionOnly as uhg_check and uhg.
    """
    uhg_check, uhg = EbeDBReader(database).getFactoryFunctions()
