from django.urls import include, path
from rest_framework import routers

from . import views

router = routers.DefaultRouter()
router.register(r'sets/csv', views.CSVSetViewSet, basename='csvSet')

urlpatterns = [
    path('', include(router.urls)),
]
