#!/usr/bin/env python
"""
    Collection of utilities related to lists.
"""

def isIterable(quantity):
    """
        Return whether "quantity" is iterable. String is not considered as
        iterable.
    """
    if type(quantity) == type(""):
        return False # string is not considered as "iterable" here
    try:
        iter(quantity)
        return True
    except TypeError:
        return False
