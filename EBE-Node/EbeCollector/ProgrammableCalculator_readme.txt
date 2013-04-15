
==================================================
    Document for ProgrammableCalculator Module
==================================================

The ProgrammableCalculator module provides a programmable calculator, which is a caluculator that accepts a string and returns a value from it. Because of the power of Python, normal expression can be easily evaluated; however the focus of this module is the ability of being "programmable". The "programmable" ability, to be more precise, is the ability of replacing special symbols by other symbols/function calls expressed in strings, and this is supported by the ProgrammableCalculator class.

The ProgrammableCalculator class uses one dictionary to specify the replacement rule, which can be set either using constructor or by an explicit call to the setRules member function.


Such replacement is performed recursively
