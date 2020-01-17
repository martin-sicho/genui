"""
urls

Created by: Martin Sicho
On: 04-12-19, 15:01
"""

from django.urls import path, include
from rest_framework import routers

import commons.views
from compounds.models import MolSet
from . import views

# Routers provide an easy way of automatically determining the URL conf.
router = routers.DefaultRouter()
router.register(r'', views.MoleculeViewSet, basename='compound')
router.register(r'sets/all', views.MolSetViewSet, basename='molset')
router.register(r'sets/chembl', views.ChEMBLSetViewSet, basename='chemblSet')

routes = [
    path('sets/all/', views.MolSetListView.as_view())
    , path('sets/<int:pk>/tasks/all/', commons.views.ModelTasksView.as_view(model_class=MolSet))
    , path('sets/<int:pk>/tasks/started/', commons.views.ModelTasksView.as_view(started_only=True, model_class=MolSet))
    , path('sets/<int:pk>/molecules/', views.MolSetMoleculesView.as_view())
]

urlpatterns = [
    path('', include(routes)),
    path('', include(router.urls)),
]
