#!/usr/bin/env python
"""
    This module implements a programmable calculator.
"""

import re

class ProgrammableCalculator(object):
    """
        This class implements the programmable calculator, which is a string
        expression evaluator with additional support to expression replacement.
        A dictionary of rules need to be set first set, which can be done either
        using the constructor or the setRules function. To see how to define a
        rule see the docstring for the applySingleRule function.
    """
    def __index__(self, rules=None):
        """
            Sets the given "rules" if any.
        """
        self.setRules(rules)
    
    def setRules(self, rules):
        """
            Stores the given rules.
        """
        self.rules = rules

    def applySingleRule(self, item, expression):
        """
            Apply a single rule to the given expression.

            A single rule is an item in the rule dictionary. Its key can either
            be a regulation pattern or a regulation string. Its value can either
            be a string, or a string containing "%s" formatting symbol, or a
            function. A replacement is triggered when the key matches, in which
            case the replacement rule is one of the following:

            1) string without formatting symbol: it will be used to replace the
                key string.
            2) string with formatting symbol: it will be formatted against the
                groups() after match then replace the key string.
            3) function: it will be called with argument groups() after match
                then the returned value of the function which should be a string
                will be used to replace the key string.

        """
        # base case checking
        key, value = item
        matched = re.match(key, expression)
        if not matched: return expression

        # matched!
        if isinstance(value, str):
            if "%" not in key:
                # case 1)
                pass
            else:
                # case 2)
                pass
        elif hasattr(value, '__call__'):
            # case 3)
            pass



    def evaluate(self, expression):
        """
            Evaluate the string argument "expression" after recursively applied
            all the rules.
        """
        pass

if __name__ == '__main__':
    import doctest
    doctest.testfile("ProgrammableCalculator_readme.txt")
