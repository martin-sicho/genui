"""genui URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path, include, re_path

from compounds.views import ChEMBLSetViewSet, MoleculeViewSet
from . import views
from rest_framework import routers
from projects.views import GenUIProjectViewSet
from django.views.generic import TemplateView
from rest_framework.schemas import get_schema_view

# Routers provide an easy way of automatically determining the URL conf.
router = routers.DefaultRouter()
router.register(r'projects', GenUIProjectViewSet, basename='project')
router.register(r'compoundSets/chembl', ChEMBLSetViewSet, basename='chemblSet')
router.register(r'compounds', MoleculeViewSet, basename='compound')

urlpatterns = [
    path('admin/', admin.site.urls),
    path('api-auth/', include('rest_framework.urls')),
    # path('qsar/', include('qsar.urls')),
    # path('compounds/', include('compounds.urls')),
    path('api/', include(router.urls)),
    path('api/schema/', get_schema_view(
        title="GenUI API",
        description="This is the schema of the root API...",
        version="1.0.0"
    ), name='openapi-schema-genui'),
    path('api/swagger-ui/', TemplateView.as_view(
        template_name='genui/swagger-ui.html',
        extra_context={'schema_url':'openapi-schema-genui'}
    ), name='swagger-ui-genui'),
    re_path(r'^api/celery-progress/', include('celery_progress.urls')),
    # path('api/projects/', include('projects.urls')),

    # if it is not a direct request to backend, serve the frontend app
    path('', views.FrontendAppView.as_view()),
    re_path(r'^(?:.*)/?$', views.FrontendAppView.as_view())
]
