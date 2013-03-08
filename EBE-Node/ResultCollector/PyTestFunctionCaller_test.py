#!/usr/bin/env python

def test_func():
    """ Dummy. """
    #print("Test function successfully called.")
    return True


if __name__ == '__main__':
    from PyTestFunctionCaller import callTestFunctionsInModule
    callTestFunctionsInModule("PyTestFunctionCaller_test")
