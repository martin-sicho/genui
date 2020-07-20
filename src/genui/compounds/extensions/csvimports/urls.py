"""
urls

Created by: Martin Sicho
On: 7/15/20, 4:36 PM
"""

from django.urls import include, path
from rest_framework import routers

from . import views

router = routers.DefaultRouter()
router.register(r'sets/csv', views.CSVSetViewSet, basename='csvSet')

urlpatterns = [
    path('', include(router.urls)),
]
