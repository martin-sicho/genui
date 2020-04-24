"""
urls

Created by: Martin Sicho
On: 25-02-20, 16:29
"""
from django.urls import path, include
from rest_framework import routers

import commons.views
import modelling.views
from maps.models import Map
from . import views
import qsar.views

router = routers.DefaultRouter()
router.register(r'algorithms',views.MappingAlgViewSet, basename='mapping-algorithm')
router.register(r'descriptors', qsar.views.DescriptorGroupsViewSet, basename='descriptor')
router.register(r'', views.MapViewSet, basename='map')


routes = [
    path('<int:pk>/tasks/all/', commons.views.ModelTasksView.as_view(model_class=Map))
    , path('<int:pk>/tasks/started/', commons.views.ModelTasksView.as_view(started_only=True, model_class=Map))
    , path('<int:pk>/files/', modelling.views.ModelFileView.as_view(model_class=Map), name="map-files-list")
    , path('<int:pk>/points/', views.PointsView.as_view(model_class=Map), name="map-points-list")
]

urlpatterns = [
    path('', include(routes)),
    path('', include(router.urls)),
]
