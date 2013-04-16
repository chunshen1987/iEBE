#!/usr/bin/env python
"""
    This module implements a programmable calculator.
"""

from StringSubstitution import StringSubstitution

class ProgrammableCalculator(StringSubstitution):
    """
        This class implements the programmable calculator, which is a string
        expression evaluator with additional support to expression replacement.
        It subclass the StringSubstitution class and provide a single evaluate
        function that performs substitution first then evaluate the resulting
        expression.
    """

    def evaluate(self, expression, verbose=False):
        """
            Evaluate the string argument "expression" after recursively applied
            all the rules. When verbose is set to True the expression after
            substitution will be printed to the screen.
        """
        toEvaluate, scans = self.applyAllRules(expression)
        if verbose: print("Evaluating:\n" + toEvaluate)
        return eval(toEvaluate)

if __name__ == '__main__':
    import doctest
    doctest.testfile("ProgrammableCalculator_readme.txt")
