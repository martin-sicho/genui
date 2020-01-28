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
from django.conf import settings
from django.conf.urls.static import static
from django.contrib import admin
from django.urls import path, include, re_path
from drf_yasg import openapi
from drf_yasg.views import get_schema_view

from . import views
from commons.views import TaskProgressView

schema_view = get_schema_view(
   openapi.Info(
      title="GenUI API",
      default_version='v0',
      description="API to interact with the GenUI backend server.",
      # terms_of_service="https://www.google.com/policies/terms/",
      # contact=openapi.Contact(email="contact@something.local"),
      # license=openapi.License(name="BSD License"), # FIXME: needs to be changed
   ),
   public=True,
)

urlpatterns = [
    path('admin/', admin.site.urls),
    path('api/auth/', include('rest_framework.urls')),
    re_path(r'^api/celery-progress/(?P<task_id>[\w-]+)/$', TaskProgressView.as_view()),
    path('api/projects/', include('projects.urls')),
    path('api/compounds/', include('compounds.urls')),
    path('api/qsar/', include('qsar.urls')),
    path('api/generators/', include('generators.urls')),
    re_path(r'^api/schema/swagger(?P<format>\.json|\.yaml)$', schema_view.without_ui(cache_timeout=0), name='schema-json'),
    re_path(r'^api/(swagger/)?$', schema_view.with_ui('swagger', cache_timeout=0), name='schema-swagger-ui'),
    re_path(r'^api/redoc/$', schema_view.with_ui('redoc', cache_timeout=0), name='schema-redoc'),

] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT) + [
    # if it is not a direct request to backend, serve the frontend app
    path('', views.FrontendAppView.as_view()),
    re_path(r'^(?:.*)/?$', views.FrontendAppView.as_view())
]
