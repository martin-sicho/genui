import traceback

from django.conf import settings
from django.db import transaction
from django.shortcuts import render

# Create your views here.
from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema
from rest_framework import pagination, mixins, viewsets, generics, status, parsers
from rest_framework.exceptions import NotFound
from rest_framework.response import Response
from rest_framework.views import APIView

import modelling.models
import modelling.serializers
from commons.views import FilterToProjectMixIn


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

class ModelFileView(
    generics.ListCreateAPIView
):
    queryset = modelling.models.ModelFile.objects.all()
    serializer_class = modelling.serializers.ModelFileSerializer

    def create(self, request, *args, **kwargs):
        request.data["model"] = self.kwargs['pk']
        serializer = self.get_serializer_class()(data=request.data)
        if serializer.is_valid():
            with transaction.atomic():
                instance = serializer.create(serializer.validated_data)
            ret = self.serializer_class(instance).data
            return Response(ret, status=status.HTTP_201_CREATED)
        else:
            print(serializer.errors)
            print(serializer.initial_data)
            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

class ModelViewSet(FilterToProjectMixIn, viewsets.ModelViewSet):
    init_serializer_class = None
    builder_class = None
    build_task = None

    def get_serializer_class(self):
        if self.action == 'create':
            if self.init_serializer_class:
                return self.init_serializer_class
            else:
                raise Exception("No init serializer specified. This is required.")
        # elif self.action in ('update', 'partial_update',):
        #     return ChEMBLSetUpdateSerializer
        else:
            return super().get_serializer_class()

    project_id_param = openapi.Parameter('project_id', openapi.IN_QUERY, description="ID of a project to limit the list of results to.", type=openapi.TYPE_NUMBER)
    @swagger_auto_schema(
        operation_description="List all models. Supply a project ID to get only models specific to a particular project."
        # , methods=['GET']
        , manual_parameters=[project_id_param]
        # , responses={200: GenericMolSetSerializer(many=True)}
    )
    def list(self, request, *args, **kwargs):
        return super().list(request, *args, **kwargs)

    def create(self, request, *args, **kwargs):
        if not self.builder_class or not self.build_task:
            raise Exception("No model builder class or build task specified. Check the viewset class definition.")

        serializer = self.get_serializer_class()(data=request.data, builder_class=self.builder_class)
        if serializer.is_valid():
            with transaction.atomic():
                instance = serializer.create(serializer.validated_data)

            task = None
            try:
                task_id = None
                if instance.build:
                    task = instance.apply_async(self.build_task, args=[instance.pk, self.builder_class.__name__])
                    task_id = task.id
                ret = self.serializer_class(instance).data
                ret["taskID"] = task_id
                return Response(ret, status=status.HTTP_201_CREATED)
            except Exception as exp:
                traceback.print_exc()
                if task and task.id:
                    settings.CURRENT_CELERY_APP.control.revoke(task_id=task.id, terminate=True)
                instance.delete()
                return Response({"error" : repr(exp)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
        print(serializer.errors)
        print(serializer.initial_data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)