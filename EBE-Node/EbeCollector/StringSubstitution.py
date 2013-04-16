#!/usr/bin/env python
"""
    This module implements a StringSubstitution class.
"""

import re

class StringSubstitution(object):
    """
        This class supplies methods for performing generalized string
        substitution using rules.
    """
    
    @staticmethod
    def applySingleRule(rule, expression, restrictNumberOfRecursionsTo=None):
        """
            Apply a single rule to the given expression.

            A single rule is an (key, value) tuple. The key can be a string or a
            compiled regulation pattern. The value can either be a string or a
            function. A substitution is triggered when the key matches, in which
            case the substitution rule is one of the following:

            1) value is a string: it will be formatted (format function) against
                the matched groups() then used for substitution.

            2) value is a function: it will be called with the matched object
                then the returned value of the function which should be a string
                will be used to perform the substitution.

            3) value is neither types: they will be converted to strings using
                str().

            If restrictNumberOfSubs is set, the substitution will be stopped
            after the specified number of recursions of substitutions are
            performed. Here one recursion means that to replace all the
            occurence (possibly many) of the first matched result (single) by
            the string generated using "value".

            This function returns the string after substitution and the number
            of recursion the substitutions are performed.

        """
        # base case checking
        key, value = rule

        # unify to functional call with matched object
        if isinstance(value, str):
            # string
            formulate = lambda matchedResult: value.format(matchedResult.groups())
        elif hasattr(value, '__call__'):
            # function
            formulate = value
        else:
            # brute force conversion to string
            formulate = lambda matchedResult: str(value)

        # loop until no more matches
        substitutedString = expression
        numberOfRecursions = 0
        matched = re.search(key, substitutedString)
        while matched:
            # generate substitution string and perform substitution
            toBeSubstituted = formulate(matched)
            substitutedString = substitutedString.replace(matched.group(), toBeSubstituted)
            # update number of recursions
            numberOfRecursions += 1
            if restrictNumberOfRecursionsTo and numberOfRecursions>=restrictNumberOfRecursionsTo: break
            # for the next iteration
            matched = re.search(key, substitutedString)

        return (substitutedString, numberOfRecursions)

    def __init__(self, rules=[]):
        """
            Sets the given "rules" if any.
        """
        self.setRules(rules)

    def setRules(self, rules):
        """
            Stores the given rules. Rules should be a list of (key, value)
            tuples. The rules are applied according to their order in the list.
            If noCompile is set to True, the rules are not compiled.
        """
        # compile and store all rules
        self.rules = []
        for key, value in rules:
            self.rules.append((re.compile(key), value))

    def applyAllRules(self, expression, restrictNumberOfRecursionsPerScan=None, restrictNumberOfScans=None):
        """
            Apply all the substitution rules to the expression. The substitution
            rules are set using either the constructor or the setRules function.

            Each rule is applied using the applySingleRule function, and
            restrictNumberOfRecursionsPerScan argument is passed to it as the
            restrictNumberOfRecursionsTo argument.

            All rules are applied in order, and repeatedly until no more changes
            occur or until restrictNumberOfScans number of repeatness reached.

            Return the expression after substitution and the number of scans
            performed. Changes are made if the returned number of scans is
            non-zero.

            If trackExpression is set to True, the function prints out all the
            expression after each substitution.
        """
        # repeatly apply all the rules
        numberOfScans = 0
        newExpression = expression
        changeFlag = True
        while changeFlag:
            # loop over all the rules
            changeFlag = False
            for rule in self.rules:
                newExpression, recursionNumber = self.applySingleRule(rule, newExpression, restrictNumberOfRecursionsTo=restrictNumberOfRecursionsPerScan)
                if recursionNumber>0: changeFlag = True
            if changeFlag: numberOfScans += 1
            # check number of scans
            if restrictNumberOfScans and numberOfScans>=restrictNumberOfScans: break

        return (newExpression, numberOfScans)



if __name__ == '__main__':
    import doctest
    doctest.testfile("StringSubstitution_readme.txt")
