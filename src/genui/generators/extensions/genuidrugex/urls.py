"""
urls

Created by: Martin Sicho
On: 5/3/20, 6:51 PM
"""
from django.urls import path, include
from rest_framework import routers

from genui.utils.extensions.tasks.views import ModelTasksView
from genui.models.views import ModelFileView, ModelPerformanceListView
from . import models
from . import views

router = routers.DefaultRouter()
router.register(r'drugex/networks', views.DrugExNetViewSet, basename='drugex_net')
router.register(r'drugex/agents', views.DrugExAgentViewSet, basename='drugex_agent')
router.register(r'drugex/environments', views.EnvironmentViewSet, basename='drugex_env')

# scoring methods
router.register(r'drugex/scorers/methods/all', views.ScoringMethodViewSet, basename='drugex_scoremethods_all')
router.register(r'drugex/scorers/methods/genuimodels', views.QSARScorerViewSet, basename='drugex_scoremethods_genuimodels')
router.register(r'drugex/scorers/methods/properties', views.PropertyScorerViewSet, basename='drugex_scoremethods_properties')

# modifiers
router.register(r'drugex/scorers/modifiers/all', views.ModifierViewSet, basename='drugex_modifiers_all')
router.register(r'drugex/scorers/modifiers/clipped', views.ClippedViewSet, basename='drugex_modifiers_clipped')
router.register(r'drugex/scorers/modifiers/hump', views.HumpViewSet, basename='drugex_modifiers_hump')

# scorers
router.register(r'drugex/scorers', views.ScorerViewSet, basename='drugex_scorers')

# generators
router.register(r'drugex/generators', views.GeneratorViewSet, basename='drugex_generators')

routes = [
    # networks
    path('drugex/networks/<int:pk>/tasks/all/', ModelTasksView.as_view(model_class=models.DrugExNet))
    , path('drugex/networks/<int:pk>/tasks/started/', ModelTasksView.as_view(started_only=True, model_class=models.DrugExNet))
    , path('drugex/networks/<int:pk>/performance/', ModelPerformanceListView.as_view(), name="drugex_net_perf_view")
    , path('drugex/networks/<int:pk>/files/', ModelFileView.as_view(model_class=models.DrugExNet), name="drugex-net-model-files-list")
] + [
    # agents
    path('drugex/agents/<int:pk>/tasks/all/', ModelTasksView.as_view(model_class=models.DrugExAgent))
    , path('drugex/agents/<int:pk>/tasks/started/', ModelTasksView.as_view(started_only=True, model_class=models.DrugExAgent))
    , path('drugex/agents/<int:pk>/performance/', ModelPerformanceListView.as_view(), name="drugex_agent_perf_view")
    , path('drugex/agents/<int:pk>/files/', ModelFileView.as_view(model_class=models.DrugExAgent), name="drugex-agent-model-files-list")
]

urlpatterns = [
    path('', include(routes)),
    path('', include(router.urls)),
]