from django.http import HttpResponse
from rest_framework import viewsets

from genui.commons.views import FilterToUserMixIn
from .serializers import ProjectSerializer
from .models import Project

def index(req):
    return HttpResponse("This is the main page of projects")

# ViewSets define the view behavior.
class GenUIProjectViewSet(FilterToUserMixIn, viewsets.ModelViewSet):
    queryset = Project.objects.all()
    serializer_class = ProjectSerializer
