import logging
from numpy import ndarray
import os

from django.http import HttpResponse
from django.template import RequestContext, loader

import bridge


def query(request):
    output_format = request.GET.get(bridge.FORMAT, "table")
    database_name = request.GET.get(bridge.DATABASE_PARAM, "")
    expression = request.GET.get(bridge.EXPR_PARAM, "")

    logging.info("Query with expr=%s, database=%s, fmt=%s",
                 expression, database_name, output_format)

    try:
        if not expression:
            raise ValueError()

        query_bridge = bridge.QueryBridge()
        query_bridge.set_database(database_name)
        results = query_bridge.evaluate_expression(expression)

        # Make results iterable; special case when it is a string.
        if type(results) == str:
            raise ValueError()
        elif type(results) != ndarray:
            results = [results]

        template = loader.get_template("query_%s.tmpl" % output_format)
        context = RequestContext(request, {
            "query_content": expression,
            "query_result": results})
        return HttpResponse(template.render(context))

    except ValueError:
        template = loader.get_template("query_error.tmpl")
        context = RequestContext(request, {
            "query_content": expression,
            "database_param": bridge.DATABASE_PARAM,
            "expr_param": bridge.EXPR_PARAM,
            "fmt_plain": "plain",
            "fmt_table": "table"})
        return HttpResponse(template.render(context))


def home(request):
    print bridge.get_available_database_names()
    template = loader.get_template("query_homepage.tmpl")
    context = RequestContext(request, {
        "database_names": bridge.get_available_database_names(),
        "formats": ["plain", "table"]})
    return HttpResponse(template.render(context))
