import traceback

from django.conf import settings
from django.db import transaction

# Create your views here.
from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema
from rest_framework import mixins, viewsets, generics, status
from rest_framework.exceptions import NotFound
from rest_framework.parsers import MultiPartParser
from rest_framework.response import Response

from genui import celery_app
from genui.utils.inspection import getFullName
from genui.utils.extensions.tasks.utils import runTask
from genui.accounts.serializers import FilterToUserMixIn
from genui.projects.serializers import FilterToProjectMixIn
from genui.utils.pagination import GenuiPagination
from genui.models.models import ModelFile, ModelPerformance, Algorithm, ModelPerformanceMetric, Model
from genui.models.serializers import ModelFileSerializer, ModelPerformanceSerializer, AlgorithmSerializer, \
    ModelPerformanceMetricSerializer


class PerformancePagination(GenuiPagination):
    page_size = 10

class FilterToModelMixin:
    lookup_field = "model"
    model_class = None

    def get_queryset(self):
        queryset = super().get_queryset()
        if "pk" in self.kwargs:
            pk = self.kwargs["pk"]
            model_class = self.model_class if self.model_class else Model
            try:
                if issubclass(self.__class__, FilterToUserMixIn) and self.request.user and not self.request.user.is_anonymous:
                    model_class.objects.get(pk=pk, project__owner=self.request.user)
                else:
                    model_class.objects.get(pk=pk)
            except model_class.DoesNotExist:
                raise NotFound(f"No model found: {pk}.", status.HTTP_400_BAD_REQUEST)
            lookup = self.lookup_field + "__id"
            return queryset.filter(**{ lookup: pk})
        else:
            return queryset

class MetricsViewSet(
    mixins.ListModelMixin,
    mixins.RetrieveModelMixin,
    viewsets.GenericViewSet
):
    queryset = ModelPerformanceMetric.objects.all()
    serializer_class = ModelPerformanceMetricSerializer


class AlgorithmViewSet(
    mixins.ListModelMixin,
    mixins.RetrieveModelMixin,
    viewsets.GenericViewSet
):
    queryset = Algorithm.objects.all()
    serializer_class = AlgorithmSerializer


class ModelPerformanceListView(
    FilterToModelMixin,
    FilterToUserMixIn,
    generics.ListAPIView
):
    queryset = ModelPerformance.objects.order_by('id')
    serializer_class = ModelPerformanceSerializer
    pagination_class = PerformancePagination
    owner_relation = "model__project__owner"

class ModelFileView(
    FilterToModelMixin,
    FilterToUserMixIn,
    generics.ListCreateAPIView
):
    parser_classes = (MultiPartParser,)
    queryset = ModelFile.objects.all()
    serializer_class = ModelFileSerializer
    lookup_field = "modelInstance"
    owner_relation = "modelInstance__project__owner"

    def create(self, request, *args, **kwargs):
        try:
            Model.objects.get(pk=self.kwargs['pk'], project__owner=request.user)
        except Model.DoesNotExist:
            return Response({"error" : f"Model does not exist: {self.kwargs['pk']}"}, status=status.HTTP_404_NOT_FOUND)

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

class PredictMixIn:
    predict_task = None

    def get_predict_task(self):
        if not self.predict_task:
            # FIXME: this should be checked in metaclass
            raise Exception("Predict task needs to be set.")
        return self.predict_task

class BuildMixIn:
    init_serializer_class = None
    builder_class = None
    build_task = None

    def get_serializer_class(self):
        if self.action == 'create':
            if self.init_serializer_class:
                return self.init_serializer_class
            else:
                raise Exception("No init serializer specified. This is required.")
        else:
            return super().get_serializer_class()

    def get_builder_class(self, asName=True):
        if not self.builder_class:
            # FIXME: this should be checked in metaclass
            raise Exception("Builder class needs to be set.")
        if asName:
            return getFullName(self.builder_class)
        else:
            return self.builder_class

    def get_build_task(self):
        if not self.build_task:
            # FIXME: this should be checked in metaclass
            raise Exception("Build task needs to be set.")
        return self.build_task

class ModelViewSet(
    FilterToProjectMixIn
    , FilterToUserMixIn
    , BuildMixIn
    , viewsets.ModelViewSet
):
    owner_relation = "project__owner"

    project_id_param = openapi.Parameter('project_id', openapi.IN_QUERY, description="ID of a project to limit the list of results to.", type=openapi.TYPE_NUMBER)
    @swagger_auto_schema(
        operation_description="List all models. Supply a project ID to get only models specific to a particular project."
        # , methods=['GET']
        , manual_parameters=[project_id_param]
        # , responses={200: GenericMolSetSerializer(many=True)}
    )
    def list(self, request, *args, **kwargs):
        return super().list(request, *args, **kwargs)

    def get_builder_kwargs(self):
        return dict()

    def create(self, request, *args, **kwargs):
        serializer = self.get_serializer_class()(data=request.data, builder_class=self.get_builder_class(asName=False))
        if serializer.is_valid():
            with transaction.atomic():
                instance = serializer.create(serializer.validated_data)

            task = None
            try:
                task_id = None
                if not hasattr(instance, "build") or instance.build:
                    task, task_id = runTask(
                        self.get_build_task(),
                        instance=instance,
                        eager=hasattr(settings, 'CELERY_TASK_ALWAYS_EAGER') and settings.CELERY_TASK_ALWAYS_EAGER,
                        args=(
                            instance.pk,
                            self.get_builder_class()
                        ),
                        kwargs=self.get_builder_kwargs()
                    )
                ret = self.serializer_class(instance).data
                ret["taskID"] = task_id
                return Response(ret, status=status.HTTP_201_CREATED)
            except Exception as exp:
                traceback.print_exc()
                if task and task.id:
                    celery_app.control.revoke(task_id=task.id, terminate=True)
                instance.delete()
                return Response({"error" : repr(exp)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
        print(serializer.errors)
        print(serializer.initial_data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)