
==================================================
    Document for ProgrammableCalculator Module
==================================================

The ProgrammableCalculator module provides a programmable calculator, which is a caluculator that accepts a string and returns a value from it. Because of the power of Python, normal expression can be easily evaluated; however the focus of this module is the ability of being "programmable". The "programmable" ability, to be more precise, is the ability of replacing special symbols by other symbols/function calls expressed in strings, and this is supported by the StringSubstitution class.

The ProgrammableCalculator(expression) class is a subclass of the StringSubstitution class and it provides one more function that evaluates an expression. As an example, the following is a calculator that evaluates an expression after replacing all "a" to 1 and all "b" to 2:
>>> from ProgrammableCalculator import ProgrammableCalculator
>>> calculator = ProgrammableCalculator( (("a",1),("b",2)) )
>>> calculator.evaluate( "a+b*4" )
9

An additional verbose argument can be supplied to print out the expression before evaluation:
>>> calculator.evaluate( "a+b*4", verbose=True)
Evaluating:
1+2*4
9

The power of the calculator resides in the design of the substitution rules, for more details see the StringSubstitution class documentation.

===========
The END
===========
