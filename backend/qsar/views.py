import traceback

from django.conf import settings
from django.db import transaction
from django.http import HttpResponse
from django.shortcuts import render
from rest_framework import viewsets, mixins, pagination, status
from rest_framework.response import Response

from qsar.tasks import buildModel
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

            task = None
            try:
                task = instance.apply_async(buildModel, args=[instance.pk, 'BasicQSARModelBuilder'])
                ret = serializers.QSARModelSerializer(instance).data
                ret["taskID"] = task.id
                return Response(ret, status=status.HTTP_201_CREATED)
            except Exception as exp:
                traceback.print_exc()
                if task and task.id:
                    settings.CURRENT_CELERY_INSTANCE.control.revoke(task_id=task.id, terminate=True)
                instance.delete()
                return Response({"error" : repr(exp)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
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
