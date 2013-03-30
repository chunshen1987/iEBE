#!/usr/bin/env python
"""
    Collection of utilities related to lists.
"""

def isIterable(quantity):
    """
        Return whether "quantity" is iterable. String is not considered as
        iterable. Examples:
        >>> isIterable([1,2])
        True
        >>> isIterable("aaa")
        False
        >>> isIterable({1:2, 3:4})
        True
    """
    if type(quantity) == type(""):
        return False # string is not considered as "iterable" here
    try:
        iter(quantity)
        return True
    except TypeError:
        return False


def assignmentsToDict(**kargs):
    """
        Convert a list of argument assignments to a dictionary. The keys are the
        strings of the arguments and the values are the assigned values. For
        example:
        >>> assignmentsToDict(x=1, y=2, z=3) == {'x':1, 'y':2, 'z':3}
        True
    """
    return kargs

def stringAssignmentsToDict(assignments):
    """
        Given a single or a list of strings, each string containing assignments
        expressions, convert all the assignements to a value dictionary which
        keys are the LHS of the assignments and values are the RHS of the
        assignments. This function can be conveniently used to phrase
        assignments expressions from command line.

        It works if the argument is a single string:
        >>> stringAssignmentsToDict("x=1 y=2 z=3") == {'x':1, 'y':2, 'z':3}
        True

        If works if the argument is a list of strings, and each string may
        contain multiple assignments:
        >>> stringAssignmentsToDict(["x=1 y=2", "z=3"]) == {'x':1, 'y':2, 'z':3}
        True
    """
    # check if assignments is a single string
    if not isIterable(assignments): assignments = [assignments]
    # join the arguments then use the assignmentsToDict function to convert
    return eval( "assignmentsToDict(%s)" % ",".join(" ".join(assignments).split()) )


if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=2)
