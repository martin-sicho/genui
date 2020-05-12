from django.http import HttpResponse
from rest_framework import viewsets

from genui.accounts.serializers import FilterToUserMixIn
from .serializers import ProjectSerializer
from .models import Project

def index(req):
    return HttpResponse("This is the main page of projects")

# ViewSets define the view behavior.
class GenUIProjectViewSet(FilterToUserMixIn, viewsets.ModelViewSet):
    queryset = Project.objects.order_by('-updated')
    serializer_class = ProjectSerializer
