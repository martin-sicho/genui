"""
urls.py

Created by: Martin Sicho
On: 04-12-19, 15:06
"""

from django.urls import path, include
from rest_framework import routers

from . import views

# Routers provide an easy way of automatically determining the URL conf.
router = routers.DefaultRouter()
router.register(r'', views.GenUIProjectViewSet, basename='project')

urlpatterns = [
    # path('', views.index, name='projects-index'),
    path('', include(router.urls)),
]