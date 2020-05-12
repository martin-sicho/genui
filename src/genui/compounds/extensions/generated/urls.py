"""
urls

Created by: Martin Sicho
On: 5/12/20, 9:52 AM
"""
from django.urls import path, include
from rest_framework import routers

from . import views

router = routers.DefaultRouter()
router.register(r'sets/generated', views.GeneratedSetViewSet, basename='generatedSet')

urlpatterns = [
    path('', include(router.urls)),
]

