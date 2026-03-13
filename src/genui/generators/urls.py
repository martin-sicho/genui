
from django.urls import path, include
from rest_framework import routers

from . import views
from genui.utils.inspection import discover_extensions_urlpatterns
from .apps import GeneratorsConfig

router = routers.DefaultRouter()
router.register(r'all', views.GeneratorViewSet, basename='generator')
router.register(r'algorithms', views.GeneratorAlgorithmViewSet, basename='generator_algorithm')

urlpatterns = [
    path('', include(router.urls)),
] + discover_extensions_urlpatterns(GeneratorsConfig.name)
