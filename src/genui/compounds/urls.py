"""
urls

Created by: Martin Sicho
On: 04-12-19, 15:01
"""

from django.urls import path, include
from rest_framework import routers

from genui.utils.extensions.tasks.views import ModelTasksView
from genui.utils.inspection import discover_extensions_urlpatterns
from . import views
from .models import MolSet
from .apps import CompoundsConfig

# Routers provide an easy way of automatically determining the URL conf.

router = routers.DefaultRouter()
router.register(r'sets/all', views.MolSetViewSet, basename='molset')
router.register(r'activity/sets', views.ActivitySetViewSet, basename='activitySet')
router.register(r'', views.MoleculeViewSet, basename='compound')

routes = [
    path('sets/<int:pk>/tasks/all/', ModelTasksView.as_view(model_class=MolSet))
    , path('sets/<int:pk>/tasks/started/', ModelTasksView.as_view(started_only=True, model_class=MolSet))
    , path('sets/<int:pk>/molecules/', views.MolSetMoleculesView.as_view(), name='moleculesInSet')
]

urlpatterns = [
    path('', include(routes)),
    path('', include(router.urls)),
] + discover_extensions_urlpatterns(CompoundsConfig.name)
