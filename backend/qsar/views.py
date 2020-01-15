from django.db import transaction
from django.http import HttpResponse
from django.shortcuts import render
from rest_framework import viewsets, mixins, pagination, status
from rest_framework.response import Response

from . import models
from . import serializers

class PerformancePagination(pagination.PageNumberPagination):
    page_size = 10

class QSARModelViewSet(viewsets.ModelViewSet):
    queryset = models.QSARModel.objects.all()
    serializer_class = serializers.QSARModelSerializer

    def get_serializer_class(self):
        if self.action == 'create':
            return serializers.QSARModelSerializerInit
        # elif self.action in ('update', 'partial_update',):
        #     return ChEMBLSetUpdateSerializer
        else:
            return super().get_serializer_class()

    def create(self, request, *args, **kwargs):
        serializer = serializers.QSARModelSerializerInit(data=request.data)
        if serializer.is_valid():
            with transaction.atomic():
                instance = serializer.create(serializer.validated_data)

            # TODO: create a task to train this model

            return Response({"status", "ok"}, status=status.HTTP_201_CREATED)
        print(serializer.errors)
        print(serializer.initial_data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

class MetricsViewSet(
    mixins.ListModelMixin,
    mixins.RetrieveModelMixin,
    viewsets.GenericViewSet
):
    queryset = models.ModelPerformanceMetric.objects.all()
    serializer_class = serializers.ModelPerformanceMetricSerializer

class AlgorithmViewSet(
    mixins.ListModelMixin,
    mixins.RetrieveModelMixin,
    viewsets.GenericViewSet
):
    queryset = models.Algorithm.objects.all()
    serializer_class = serializers.AlgorithmSerializer

class DescriptorGroupsViewSet(
    mixins.ListModelMixin,
    mixins.RetrieveModelMixin,
    viewsets.GenericViewSet
):
    queryset = models.DescriptorGroup.objects.all()
    serializer_class = serializers.DescriptorGroupSerializer

class ModelPerformanceViewSet(
    mixins.ListModelMixin,
    viewsets.GenericViewSet
):
    queryset = models.ModelPerformance.objects.all()
    serializer_class = serializers.ModelPerformanceSerializer
    pagination_class = PerformancePagination
