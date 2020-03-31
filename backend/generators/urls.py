"""
urls

Created by: Martin Sicho
On: 02-12-19, 17:18
"""

from django.urls import path, include
from rest_framework import routers

import commons.views
import modelling.views
from .models import DrugExNet, DrugExAgent
from . import views

router = routers.DefaultRouter()
router.register(r'all', views.GeneratorViewSet, basename='generator')
router.register(r'drugex/networks', views.DrugExNetViewSet, basename='drugex_net')
router.register(r'drugex/agents', views.DrugExAgentViewSet, basename='drugex_agent')
router.register(r'algorithms', views.GeneratorAlgorithmViewSet, basename='generator_algorithm')
router.register(r'metrics', views.GeneratorMetricsViewSet, basename='generator_metric')


routes = [
    path('drugex/networks/<int:pk>/tasks/all/', commons.views.ModelTasksView.as_view(model_class=DrugExNet))
    , path('drugex/networks/<int:pk>/tasks/started/', commons.views.ModelTasksView.as_view(started_only=True, model_class=DrugExNet))
    , path('drugex/networks/<int:pk>/performance/', modelling.views.ModelPerformanceListView.as_view(), name="drugex_net_perf_view")
    , path('drugex/networks/<int:pk>/files/', modelling.views.ModelFileView.as_view(model_class=DrugExNet), name="drugex-net-model-files-list")
] + [
    path('drugex/agents/<int:pk>/tasks/all/', commons.views.ModelTasksView.as_view(model_class=DrugExAgent))
    , path('drugex/agents/<int:pk>/tasks/started/', commons.views.ModelTasksView.as_view(started_only=True, model_class=DrugExAgent))
    , path('drugex/agents/<int:pk>/performance/', modelling.views.ModelPerformanceListView.as_view(), name="drugex_agent_perf_view")
    , path('drugex/agents/<int:pk>/files/', modelling.views.ModelFileView.as_view(model_class=DrugExAgent), name="drugex-agent-model-files-list")
]

urlpatterns = [
    path('', include(routes)),
    path('', include(router.urls)),
]
