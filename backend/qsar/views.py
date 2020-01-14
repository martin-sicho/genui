from django.http import HttpResponse
from django.shortcuts import render
from rest_framework import viewsets, mixins, pagination

from commons.helpers import getSubclassesFromModule
from . import models
from . import serializers

from .algorithms import bases
from .algorithms import algorithms as algs
from .algorithms import metrics as metrics

class PerformancePagination(pagination.PageNumberPagination):
    page_size = 10

class QSARModelViewSet(viewsets.ModelViewSet):
    queryset = models.QSARModel.objects.all()
    serializer_class = serializers.QSARModelSerializer

class MetricsViewSet(
    mixins.ListModelMixin,
    mixins.RetrieveModelMixin,
    viewsets.GenericViewSet
):
    METRICS_ALL = [x.getDjangoModel() for x in getSubclassesFromModule(bases.ValidationMetric, metrics)]

    queryset = models.ModelPerformanceMetric.objects.all()
    serializer_class = serializers.ModelPerformanceMetricSerializer

class AlgorithmViewSet(
    mixins.ListModelMixin,
    mixins.RetrieveModelMixin,
    viewsets.GenericViewSet
):
    ALGORITHMS_ALL = [x.getDjangoModel() for x in getSubclassesFromModule(bases.BaseAlgorithm, algs)]

    queryset = models.Algorithm.objects.all()
    serializer_class = serializers.AlgorithmSerializer

class ModelPerformanceViewSet(
    mixins.ListModelMixin,
    viewsets.GenericViewSet
):
    queryset = models.ModelPerformance.objects.all()
    serializer_class = serializers.ModelPerformanceSerializer
    pagination_class = PerformancePagination
