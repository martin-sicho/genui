
from django.urls import path, include
from rest_framework import routers

from genui.utils.extensions.tasks.views import ModelTasksView
from genui.models.views import ModelFileView, ModelPerformanceListView
from .models import QSARModel
from . import views
from genui.utils.inspection import discover_extensions_urlpatterns
from .apps import QsarConfig

router = routers.DefaultRouter()
router.register(r'models', views.QSARModelViewSet, basename='model')
router.register(r'algorithms',views.QSARAlgorithmViewSet, basename='algorithm')
router.register(r'embeddings', views.EmbeddingGroupsViewSet, basename='embedding')
router.register(r'hyper-parameters', views.QSARHyperParameterOptimizationViewSet, basename='hyper-parameter')
# router.register(r'predictions', views.ModelPredictionsViewSet, basename='prediction')


routes = [
    path('models/<int:pk>/tasks/all/', ModelTasksView.as_view(model_class=QSARModel))
    , path('models/<int:pk>/tasks/started/', ModelTasksView.as_view(started_only=True, model_class=QSARModel))
    , path('models/<int:pk>/performance/', ModelPerformanceListView.as_view())
    , path('models/<int:pk>/files/', ModelFileView.as_view(model_class=QSARModel), name="qsar-model-files-list")
    , path('models/qsprpred/sklearn/', views.QSPRPredSklearnModelViewSet.as_view(({'get': 'list'})))
    , path('models/qsprpred/sklearn/mode/<str:mode>/', views.QSPRPredSklearnModelViewSet.as_view(({'get': 'mode_model'})))
    , path('models/qsprpred/sklearn/<str:algorithm>/', views.QSPRPredSklearnModelViewSet.as_view(({'get': 'retrieve'})))
    , path('models/qsprpred/sklearn/<str:algorithm>/type/', views.QSPRPredSklearnModelViewSet.as_view(({'get': 'get_type'})))
    , path('models/qsprpred/sklearn/<str:algorithm>/params/', views.QSPRPredSklearnModelViewSet.as_view(({'get': 'get_params'})))
    , path('metrics/<str:mode>/', views.list_metrics_by_mode)
    , path('metrics/', views.list_metrics)
    , path('data-splits/', views.list_data_splits)
    , path('data-splits/<str:data_split>/params', views.get_data_splits_params)
    , path('aggregation-functions/', views.list_aggregation_functions)
    , path('data-splits/scaffolds/', views.list_scaffolds)
    , path('data-splits/clustering/', views.list_clustering)
]

urlpatterns = [
    path('', include(routes)),
    path('', include(router.urls)),
] + discover_extensions_urlpatterns(QsarConfig.name)
