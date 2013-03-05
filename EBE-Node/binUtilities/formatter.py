#!/usr/bin/env python
""" All functions in this module are related to data and "formats".
    A "format" is a dictionary that defines an order. For example,
    the format {"a":1, "b":2} specifying that data of type "a"
    should be in the 1st column and data of type "b" should be in
    the 2nd column. The order of the items in the dictionary is
    meaningless: {"b":2, "a":1} specifies the same order. The column
    index starts with 0.
"""

# Version 1

def getReorderingFunction(from_format, to_format, indexStartsWith1=True):
    """ Both from_format and to_format are "format" dictionaries.
    This function returns a function that, when acting on a list
    with format "from_format", transforms it into a list of format
    "to_format". For example, if from_format = {"a":1, "b":2} and
    to_format = {"b":1} then the returned function will be such that
    it returns a list whose 1st element is the 2nd element of its
    input. When indexStartsWith1 is set to True all indices are assumed
    to start with 1; otherwise 0.
    """
    if indexStartsWith1:
        shift = -1 # -1 to convert to indices that start with 0
    else:
        shift = 0

    # first sort keys in to_format to the correct order:
    to_format_items = to_format.items()
    to_format_items.sort(key=lambda var:var[1])
    to_format_sorted_keys = [var[0] for var in to_format_items]

    # next get the desired column index list from from_format
    desired_column_order = [int(from_format[var])+shift for var in to_format_sorted_keys]

    # finally build and return a function that does the reordering
    return lambda var_input: [var_input[var] for var in desired_column_order]
