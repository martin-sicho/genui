"""
urls

Created by: Martin Sicho
On: 02-12-19, 17:18
"""

from django.urls import path, include
from rest_framework import routers

from genui.extensions.tasks.views import ModelTasksView
from genui.modelling.views import ModelFileView, ModelPerformanceListView
from .models import QSARModel
from . import views
from genui.extensions.utils import discover_extensions_urlpatterns
from .apps import QsarConfig

router = routers.DefaultRouter()
router.register(r'models', views.QSARModelViewSet, basename='model')
router.register(r'algorithms',views.QSARAlgorithmViewSet, basename='algorithm')
router.register(r'metrics', views.QSARMetricsViewSet, basename='metric')
router.register(r'descriptors', views.DescriptorGroupsViewSet, basename='descriptor')
# router.register(r'predictions', views.ModelPredictionsViewSet, basename='prediction')


routes = [
    path('models/<int:pk>/tasks/all/', ModelTasksView.as_view(model_class=QSARModel))
    , path('models/<int:pk>/tasks/started/', ModelTasksView.as_view(started_only=True, model_class=QSARModel))
    , path('models/<int:pk>/performance/', ModelPerformanceListView.as_view())
    , path('models/<int:pk>/files/', ModelFileView.as_view(model_class=QSARModel), name="qsar-model-files-list")
]

urlpatterns = [
    path('', include(routes)),
    path('', include(router.urls)),
] + discover_extensions_urlpatterns(QsarConfig.name)
