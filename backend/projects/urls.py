"""
urls.py

Created by: Martin Sicho
On: 04-12-19, 15:06
"""

from django.urls import path, include
from django.views.generic import TemplateView
from rest_framework import routers
from rest_framework.schemas import get_schema_view

from . import views

# Routers provide an easy way of automatically determining the URL conf.
router = routers.DefaultRouter()
router.register(r'', views.GenUIProjectViewSet, basename='project')

urlpatterns = [
    # path('', views.index, name='projects-index'),
    path('', include(router.urls)),
    path('schema', get_schema_view(
        title="Your Project",
        description="API for all things …",
        version="1.0.0"
    ), name='openapi-schema'),
    path('swagger-ui/', TemplateView.as_view(
        template_name='projects/swagger-ui.html',
        extra_context={'schema_url':'openapi-schema'}
    ), name='swagger-ui')
]