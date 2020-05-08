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
import urllib.parse
from django.conf import settings
from django.conf.urls.static import static
from django.urls import re_path, path
from django.views.generic import RedirectView
from drf_yasg import openapi
from drf_yasg.views import get_schema_view

from genui.utils.inspection import discover_apps_urls, discover_extensions_urlpatterns
from genui.apps import API_APPS, BASE_APPS

schema_view = get_schema_view(
    openapi.Info(
        title="GenUI API",
        default_version='v0',
        description="API to interact with the GenUI backend server.",
        # terms_of_service="https://www.google.com/policies/terms/",
        # contact=openapi.Contact(email="contact@something.local"),
        # license=openapi.License(name="BSD License"), # FIXME: needs to be changed
    ),
    public=False,
    url=settings.GENUI_SETTINGS['HOST_URL'] if 'HOST_URL' in settings.GENUI_SETTINGS else None
)

base_urls = discover_apps_urls(BASE_APPS)
api_urls = discover_apps_urls(API_APPS, prefix='api', app_names_as_root=True)
api_urls += [
    re_path(r'^api/schema/swagger(?P<format>\.json|\.yaml)$', schema_view.without_ui(cache_timeout=0), name='schema-json'),
    re_path(r'^api/redoc/$', schema_view.with_ui('redoc', cache_timeout=0), name='schema-redoc'),
    re_path(r'^api/(swagger/)?$', schema_view.with_ui('swagger', cache_timeout=0), name='schema-swagger-ui'),
]
extensions_urls = discover_extensions_urlpatterns('genui')

urlpatterns = base_urls + api_urls + extensions_urls
if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

# redirect from root to the app if the path is available
# TODO: create some simple introduction page at the root (link to repo, docs, etc.)
if 'FRONTEND_APP_PATH' in settings.GENUI_SETTINGS and settings.GENUI_SETTINGS['FRONTEND_APP_PATH'] is not None:
    urlpatterns.append(path('', RedirectView.as_view(
        url=urllib.parse.urljoin(settings.GENUI_SETTINGS['HOST_URL'], settings.GENUI_SETTINGS['FRONTEND_APP_PATH']) if 'HOST_URL' in settings.GENUI_SETTINGS else settings.GENUI_SETTINGS['FRONTEND_APP_PATH']
    )))
