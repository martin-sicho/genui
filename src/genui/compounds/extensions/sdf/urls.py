"""
urls

Created by: Martin Sicho
On: 7/13/20, 1:08 PM
"""

from django.urls import include, path
from rest_framework import routers

from . import views

router = routers.DefaultRouter()
router.register(r'sets/sdf', views.SDFSetViewSet, basename='sdfSet')

urlpatterns = [
    path('', include(router.urls)),
]