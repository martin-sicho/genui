from django.shortcuts import render

# Create your views here.
from rest_framework import pagination, mixins, viewsets, generics, status
from rest_framework.exceptions import NotFound

import modelling.models
import modelling.serializers
from qsar import models


class PerformancePagination(pagination.PageNumberPagination):
    page_size = 10


class MetricsViewSet(
    mixins.ListModelMixin,
    mixins.RetrieveModelMixin,
    viewsets.GenericViewSet
):
    queryset = modelling.models.ModelPerformanceMetric.objects.all()
    serializer_class = modelling.serializers.ModelPerformanceMetricSerializer


class AlgorithmViewSet(
    mixins.ListModelMixin,
    mixins.RetrieveModelMixin,
    viewsets.GenericViewSet
):
    queryset = modelling.models.Algorithm.objects.all()
    serializer_class = modelling.serializers.AlgorithmSerializer


class ModelPerformanceListView(
   generics.ListAPIView
):
    queryset = modelling.models.ModelPerformance.objects.order_by('id')
    serializer_class = modelling.serializers.ModelPerformanceSerializer
    pagination_class = PerformancePagination

    def get_queryset(self):
        queryset = super().get_queryset()
        if "pk" in self.kwargs:
            pk = self.kwargs["pk"]
            try:
                modelling.models.Model.objects.get(pk=pk)
            except modelling.models.Model.DoesNotExist:
                raise NotFound(f"The mode with id={pk} does not exist.", status.HTTP_400_BAD_REQUEST)
            return queryset.filter(model__id=pk)
        else:
            return queryset