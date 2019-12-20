"""
urls

Created by: Martin Sicho
On: 04-12-19, 15:01
"""

from django.urls import path, include
from rest_framework import routers

from . import views

# Routers provide an easy way of automatically determining the URL conf.
router = routers.DefaultRouter()
router.register(r'', views.MoleculeViewSet, basename='compound')
router.register(r'sets/chembl', views.ChEMBLSetViewSet, basename='chemblSet')

routes = [
    path('sets/<int:pk>/tasks/all/', views.MolSetTasksView.as_view())
    , path('sets/<int:pk>/tasks/started/', views.MolSetTasksView.as_view(started_only=True))
]

urlpatterns = [
    path('', include(router.urls)),
    path('', include(routes)),
]
