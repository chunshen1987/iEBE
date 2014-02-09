from urlparse import parse_qs
from django.shortcuts import render
from django.http import HttpResponse


def query(request, params={}):
    print(params)
    return HttpResponse("Hello, world. You're at the poll index.")


def home(request):
    return HttpResponse("Hello, this is the homepage.")
