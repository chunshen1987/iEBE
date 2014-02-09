from django.shortcuts import render
from django.http import HttpResponse

import bridge

query_bridge = bridge.QueryBridge()

print
print "Hello"
print


def query(request):
    database_name = request.GET.get(bridge.DATABASE_PARAM, "")
    expression = request.GET.get(bridge.EXPR_PARAM, "")

    if not expression:
        return HttpResponse("Syntax: /query?%s=database_name&%s=expression" %
                            (bridge.DATABASE_PARAM, bridge.EXPR_PARAM))
    else:
        print database_name, expression
        query_bridge.set_database(database_name)
        return HttpResponse(query_bridge.evaluate_expression(expression))


def home(unused_request):
    return HttpResponse("Hello, this is the homepage.")
