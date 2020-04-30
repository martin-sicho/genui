"""
urls

Created by: Martin Sicho
On: 4/30/20, 7:59 PM
"""
from django.urls import path, include, re_path
from django.contrib import admin
from allauth.account.views import ConfirmEmailView
from django.conf import settings

urlpatterns = [
    path('admin/', admin.site.urls),
    path('accounts/', include('allauth.urls')),
    re_path(r'^accounts/registration/account-confirm-email/(?P<key>[-:\w]+)/$', ConfirmEmailView.as_view(),
    name='account_confirm_email'),
    path(f'api/{settings.REST_FRAMEWORK["URLS_ROOT"]}', include('rest_framework.urls')),
    path('api/accounts/', include('rest_auth.urls')),
    path('api/accounts/registration/', include('rest_auth.registration.urls')),
]