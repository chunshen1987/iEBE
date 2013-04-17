#!/usr/bin/env python
"""
    This module setups a environment for later convenient calculations using the
    EbeDBReader class.
    
    One way to use this module from command line to evaluate a single
    expression, for example:
    uhg.py database_filename "e_2(ed)"
    
    Another way to use this module by importing all its contents to an
    interactive shell, for example:
    python -ic "from uhg import *"

    In the interactive mode use the "use" function to connect to a database, use
    "h" function to print out a short help, and use the "e" function to evaluate
    an expression.
    
"""
from numpy import *
from EbeCollector import EbeDBReader

_storedEbeDBReader = None

e = lambda s: _storedEbeDBReader.evaluateExpressionOnly(s)

def use(database):
    """
        Create a EbeDBReader object and link the factory functions for
        evaluateExpression and evaluateExpressionOnly as uhg_check and uhg.
    """
    global _storedEbeDBReader
    _storedEbeDBReader = EbeDBReader(database)

def info():
    """
        Print out number of events information for all particles.
    """
    global _storedEbeDBReader
    print("Total number of events: {}".format(_storedEbeDBReader.getNumberOfEvents()))
    print("-"*60)
    print("\t{:<30}{:^20}".format("Particle","Number of events"))
    for aParticle, numberOfEvents in _storedEbeDBReader.getAttendance():
        if numberOfEvents>0: print("\t{:<30}{:^20}".format(aParticle, numberOfEvents))

def h():
    """
        Display a short help message.
    """
    print(r"""
----------------------------------
    Welcome to the UHG Module!
----------------------------------
1) Connect to a SQLite database using "use(database)" function.
2) Use "h()" function to print this message.
3) Use "info()" function to list particle storage info in the database.
4) Evaluate expressions using "e(expression)" function. The following are examples of some supported symbols, for a complete list see EbeCollector_readme.txt.

    Ecc_{m,n}(ed), E_n(sd), ecc_{m,n}(e), Phi_{m,n}(s), {r^m}(ed),
    V_{n}(pion), v_n(kaon), Psi_n(nucleon), dN/dy(pion)
    V_{n}(pT)(pion), v_n([pT1, pT2, ...])(kaon), Psi_n(pT)(pion), dN/(dydpT)(total),
    v_n[2](pion), e_n[4](ed), <symbol>, |symbol|
    
5) Some concrete examples:
    
    E_2(s), e_{3,1}(e), ecc_3(ed), Phi_2(e),
    V_2(pion), v_3(kaon), Psi_4(total), N(total),
    v_3(0.5)(pion), v_2(linspace(0,3,50))(proton)
    dN/dydpT([0,0.5,1])(total),
    v_2[2](pion), v_2[2](0.5)(kaon), e_2[2](e), e_3[4](s)

Enjoy!
          """)

if __name__ == '__main__':
    from sys import argv
    try:
        databaseFilename = argv[1]
        use(databaseFilename)
        expr = " ".join(argv[2:])
        if not expr: raise ValueError()
        print(e(expr))
    except:
        print("Usage: uhg.py database_filename 'symbols to be evaluated'")
else:
    h()
