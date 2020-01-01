"""
urls

Created by: Martin Sicho
On: 04-12-19, 15:01
"""

from django.urls import path, include
from drf_yasg.utils import swagger_auto_schema
from rest_framework import routers
from drf_yasg import openapi

from . import views

# Routers provide an easy way of automatically determining the URL conf.
from .serializers import GenericMolSetSerializer

router = routers.DefaultRouter()
router.register(r'', views.MoleculeViewSet, basename='compound')
router.register(r'sets/chembl', views.ChEMBLSetViewSet, basename='chemblSet')

project_id_param = openapi.Parameter('project_id', openapi.IN_QUERY, description="Return compound sets related to just this project.", type=openapi.TYPE_NUMBER)

routes = [
    path('sets/generic/', swagger_auto_schema(
        operation_description="List all compound sets. Can give a project ID to filter on."
        , methods=['GET']
        , manual_parameters=[project_id_param]
        , responses={200: GenericMolSetSerializer(many=True)}
)(views.MolSetListView.as_view()))
    , path('sets/<int:pk>/tasks/all/', views.MolSetTasksView.as_view())
    , path('sets/<int:pk>/tasks/started/', views.MolSetTasksView.as_view(started_only=True))
    , path('sets/<int:pk>/molecules/', views.MolSetMoleculesView.as_view())
]

urlpatterns = [
    path('', include(router.urls)),
    path('', include(routes)),
]
