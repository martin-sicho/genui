"""
urls

Created by: Martin Sicho
On: 4/27/20, 10:25 PM
"""
from django.urls import include, path
from rest_framework import routers

from . import views

router = routers.DefaultRouter()
router.register(r'sets/chembl', views.ChEMBLSetViewSet, basename='chemblSet')
router.register(r'sets/chembl/assays', views.ChEMBLAssayViewSet, basename='chemblSetAssay')
router.register(r'sets/chembl/targets', views.ChEMBLTargetViewSet, basename='chemblSetTarget')

urlpatterns = [
    path('', include(router.urls)),
]