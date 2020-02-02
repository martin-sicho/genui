"""
urls

Created by: Martin Sicho
On: 02-12-19, 17:18
"""

from django.urls import path, include
from rest_framework import routers

import commons.views
import modelling.views
from .models import DrugExNet
from . import views

router = routers.DefaultRouter()
router.register(r'drugex/networks', views.DrugExNetViewSet, basename='drugex_net')
router.register(r'drugex/agents', views.DrugExAgentViewSet, basename='drugex_agent')
router.register(r'algorithms', views.GeneratorAlgorithmViewSet, basename='generator_algorithm')

# TODO: in this case, not all metrics will likely be available so we should add some way for algorithms to specify the metrics they are compatible with (or maybe we should expose builder classes too and set it there?)
# router.register(r'metrics', modelling.views.MetricsViewSet, basename='generator_metric')


routes = [
    path('drugex/<int:pk>/tasks/all/', commons.views.ModelTasksView.as_view(model_class=DrugExNet))
    , path('drugex/<int:pk>/tasks/started/', commons.views.ModelTasksView.as_view(started_only=True, model_class=DrugExNet))
    , path('drugex/<int:pk>/performance/', modelling.views.ModelPerformanceListView.as_view(), name="drugex_perf_view")
]

urlpatterns = [
    path('', include(routes)),
    path('', include(router.urls)),
]
