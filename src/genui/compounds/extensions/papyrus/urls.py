from django.urls import include, path
from rest_framework import routers

from . import views

router = routers.DefaultRouter()
router.register(r'sets/papyrus', views.PapyrusSetViewSet, basename='papyrusSet')
router.register(r'sets/papyrus/assays', views.PapyrusAssayViewSet, basename='papyrusSetAssay')
router.register(r'sets/papyrus/targets', views.PapyrusTargetViewSet, basename='papyrusSetTarget')

urlpatterns = [
    path('', include(router.urls)),
]