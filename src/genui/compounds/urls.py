"""
urls

Created by: Martin Sicho
On: 04-12-19, 15:01
"""

from django.urls import path, include
from rest_framework import routers

from genui.commons.views import ModelTasksView
from genui.compounds.models import MolSet
from . import views

# Routers provide an easy way of automatically determining the URL conf.
router = routers.DefaultRouter()
router.register(r'sets/all', views.MolSetViewSet, basename='molset')
router.register(r'sets/chembl', views.ChEMBLSetViewSet, basename='chemblSet')
router.register(r'sets/chembl/assays', views.ChEMBLAssayViewSet, basename='chemblSetAssay')
router.register(r'sets/chembl/targets', views.ChEMBLTargetViewSet, basename='chemblSetTarget')
router.register(r'sets/generated', views.GeneratedSetViewSet, basename='generatedSet')
router.register(r'activity/sets', views.ActivitySetViewSet, basename='activitySet')
router.register(r'', views.MoleculeViewSet, basename='compound')

routes = [
    path('sets/<int:pk>/tasks/all/', ModelTasksView.as_view(model_class=MolSet))
    , path('sets/<int:pk>/tasks/started/', ModelTasksView.as_view(started_only=True, model_class=MolSet))
    , path('sets/<int:pk>/molecules/', views.MolSetMoleculesView.as_view(), name='moleculesInSet')
]

urlpatterns = [
    path('', include(routes)),
    path('', include(router.urls)),
]
