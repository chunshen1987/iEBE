from django.shortcuts import render
from django.http import HttpResponse


def query(request):
    print request
    return HttpResponse("Hello, world. You're at the poll index.")


def home(request):
    print request
    return HttpResponse("Hello, world. You're at the poll index.")
