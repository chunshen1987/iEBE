#!/usr/bin/env python
"""
    This module consists of functions dealing with collection and processing of
    data files generated during e-by-e calculations into databases.

    ***** eccentricity *****:
    This table has the format:
    event_id (int), harmonic_order (smallint), ecc_real (double), ecc_imag (double)

"""

from os import path
from sys import argv
from DBR import SqliteDB

import numpy as np


class EbeCollector:
    pass



#------------------------------------------------------------------------------
if __name__ == "__main__":
    print("----------------------------------------------------")
    print("Welcome! -- Zhi Qiu")
    print("----------------------------------------------------")
    if len(argv) == 1:
        print("Use one of the following:\n")
        print("See the source file.")
    else:
        print("Executing: "+argv[1]+"('"+"','".join(argv[2:])+"')\n")
        exec(argv[1]+"('"+"','".join(argv[2:])+"')")
