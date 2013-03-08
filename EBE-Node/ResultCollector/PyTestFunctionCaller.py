#!/usr/bin/env python

from sys import argv

testFunctionStartWith = "test_"

def callTestFunctionsInModule(module, interceptException=True):
    """
        Look for all functions whose names start with the string "test_" and call
        them. Test functions of such should have no input arguments and should
        return a boolean value indicating the successfullness. Print a
        comfirmation and its return value. Only test functions intended to be
        called should start with the symbol "test_". The argument "module" is
        the name of the module to be tested, not the filename.
    """
    exec("import %s" % module)
    exec("listOfSymbols = dir(%s)" % module)
    print("="*15 + " " + "Start testing module %s" % module + " " + "="*15)
    for aSymbol in listOfSymbols:
        if aSymbol.startswith(testFunctionStartWith):
            # it is a test function, call it
            try:
                exec("returnValue = %s.%s()" % (module, aSymbol))
            except AssertionError: # if assertion is wrong, catch it to continue the test
                returnValue = False
                if not interceptException:
                    raise
            # print out test results
            if returnValue:
                print(" "*5 + "Test %s passed." % aSymbol)
            else:
                print("-"*3+"> " + "Test %s failed!" % aSymbol)


if __name__ == '__main__':
    if len(argv)<2:
        print("Usage: PyTestFunctionCaller module_name [should_intercept_exception]")
    else:
        if len(argv)==3:
            callTestFunctionsInModule(argv[1], eval(argv[2]))
        else:
            callTestFunctionsInModule(argv[1])
