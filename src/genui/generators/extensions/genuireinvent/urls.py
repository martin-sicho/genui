# genui/src/genui/generators/extensions/genuireinvent/urls.py

from django.urls import path, include
from rest_framework import routers

from genui.utils.extensions.tasks.views import ModelTasksView
from genui.models.views import ModelFileView, ModelPerformanceListView

from . import models, views

router = routers.DefaultRouter()

router.register(r"reinvent/networks", views.ReinventNetViewSet, basename="reinvent-net")
router.register(r"reinvent/agents", views.ReinventAgentViewSet)
router.register(r"reinvent/environments", views.ReinventEnvironmentViewSet)

router.register(r"reinvent/diversity-filters", views.ReinventDiversityFilterViewSet)
router.register(r"reinvent/property-scorers", views.PropertyScorerViewSet, basename="reinvent-property-scorer")
router.register(r"reinvent/model-scorers", views.GenUIModelScorerViewSet, basename="reinvent-model-scorer")
router.register(r"reinvent/unwanted-smarts", views.UnwantedSmartsScorerViewSet, basename="reinvent-unwanted-smarts")
router.register(r"reinvent/agent-training", views.ReinventAgentTrainingViewSet)
router.register(r"reinvent/agent-validation", views.ReinventAgentValidationViewSet)
router.register(r"reinvent/runs", views.ReinventViewSet)
router.register(r"reinvent/stages", views.ReinventStageViewSet)
router.register(r"reinvent/performance", views.ModelPerformanceReinventViewSet, basename="reinvent-performance")
router.register(r"reinvent/model-files", views.ReinventModelFileViewSet, basename="reinvent-model-files")
router.register(r"reinvent/metrics", views.ReinventMetricsViewSet, basename="reinvent-metrics")

routes = [
    path("reinvent/networks/<int:pk>/tasks/all/",ModelTasksView.as_view(model_class=models.ReinventNet)),
    path("reinvent/networks/<int:pk>/tasks/started/",ModelTasksView.as_view(started_only=True, model_class=models.ReinventNet)),
    path("reinvent/networks/<int:pk>/performance/",ModelPerformanceListView.as_view(),name="reinvent_net_perf_view"),
    path("reinvent/networks/<int:pk>/files/",ModelFileView.as_view(model_class=models.ReinventNet),name="reinvent-net-model-files-list"),
    path("reinvent/runs/<int:pk>/tasks/all/", ModelTasksView.as_view(model_class=models.Reinvent)),
    path("reinvent/runs/<int:pk>/tasks/started/", ModelTasksView.as_view(started_only=True, model_class=models.Reinvent)),
]

urlpatterns = [
    path("", include(routes)),
    path("", include(router.urls)),
]