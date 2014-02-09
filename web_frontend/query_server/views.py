from django.shortcuts import render
from django.http import HttpResponse


def query(request):
    expression = request.GET.get("expr", "")
    if not expression:
        return HttpResponse("Syntax: /query?expr=expression")
    else:
        print expression
        return HttpResponse(expression)


def home(request):
    return HttpResponse("Hello, this is the homepage.")
