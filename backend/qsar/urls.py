"""
urls

Created by: Martin Sicho
On: 02-12-19, 17:18
"""

from django.urls import path, include
from rest_framework import routers

import commons.views
import modelling.views
from qsar.models import QSARModel
from . import views

router = routers.DefaultRouter()
router.register(r'models', views.QSARModelViewSet, basename='model')
router.register(r'algorithms', modelling.views.AlgorithmViewSet, basename='algorithm')
router.register(r'metrics', modelling.views.MetricsViewSet, basename='metric')
router.register(r'descriptors', views.DescriptorGroupsViewSet, basename='descriptor')


routes = [
    path('models/<int:pk>/tasks/all/', commons.views.ModelTasksView.as_view(model_class=QSARModel))
    , path('models/<int:pk>/tasks/started/', commons.views.ModelTasksView.as_view(started_only=True, model_class=QSARModel))
    , path('models/<int:pk>/performance/', modelling.views.ModelPerformanceListView.as_view())
]

urlpatterns = [
    path('', include(routes)),
    path('', include(router.urls)),
]
