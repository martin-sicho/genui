from django.http import HttpResponse
from django.shortcuts import render

def index(req):
    return HttpResponse("This is the main page of compounds")
