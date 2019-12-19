from django.http import HttpResponse
from django.shortcuts import render
from rest_framework import viewsets

from .serializers import GenUIProjectSerializer
from .models import Project


def index(req):
    return HttpResponse("This is the main page of projects")

# ViewSets define the view behavior.
class GenUIProjectViewSet(viewsets.ModelViewSet):
    queryset = Project.objects.all()
    serializer_class = GenUIProjectSerializer
    paginator = None
