import traceback

from django.db import transaction
from django.conf import settings
from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema
from rest_framework import viewsets, pagination, mixins, status, generics
from rest_framework.decorators import action
from rest_framework.response import Response
from rest_framework.viewsets import GenericViewSet
from rest_framework.schemas.openapi import AutoSchema

from genui.generators.models import GeneratedMolSet
from genui.generators.serializers import GeneratedSetInitSerializer, GeneratedSetSerializer
from genui.commons.views import FilterToProjectMixIn, FilterToUserMixIn
from genui.compounds.initializers.generated import GeneratedSetInitializer
from genui.compounds.serializers import MoleculeSerializer, MolSetSerializer, \
    GenericMolSetSerializer, ActivitySetSerializer, ActivitySerializer, \
    ActivitySetSummarySerializer
from .models import Molecule, MolSet, ActivitySet, Activity
from .tasks import populateMolSet, updateMolSet

from django_rdkit import models as djrdkit

from genui.commons.helpers import getFullName
from ..commons.tasks import runTask


class MoleculePagination(pagination.PageNumberPagination):
    page_size = 5

class ActivityPagination(pagination.PageNumberPagination):
    page_size = 10

class BaseMolSetViewSet(
    FilterToProjectMixIn
    , FilterToUserMixIn
    , viewsets.ModelViewSet
):
    class Schema(MolSetSerializer.AutoSchemaMixIn, AutoSchema):
        pass
    schema = Schema()
    initializer_class = None
    updater_class = None
    owner_relation = "project__owner"

    def get_initializer_class(self):
        if not self.initializer_class:
            raise Exception("Initializer class needs to be set.")
        return getFullName(self.initializer_class)

    def get_updater_class(self):
        if not self.updater_class:
            return self.get_initializer_class()
        return getFullName(self.updater_class)

    def get_initializer_additional_arguments(self, validated_data):
        return dict()

    def get_updater_additional_arguments(self, validated_data):
        return dict()

    project_id_param = openapi.Parameter('project_id', openapi.IN_QUERY, description="Return compound sets related to just the project with a given ID.", type=openapi.TYPE_NUMBER)
    @swagger_auto_schema(
        operation_description="List all compound sets. Can give a project ID to filter on."
        # , methods=['GET']
        , manual_parameters=[project_id_param]
    )
    def list(self, request, *args, **kwargs):
        return super().list(request, *args, **kwargs)

    def create(self, request, *args, **kwargs):
        serializer_class = self.get_serializer_class()
        serializer = serializer_class(data=request.data)
        if serializer.is_valid():
            with transaction.atomic():
                instance = serializer.create(serializer.validated_data)

            task = None
            try:
                task, task_id = runTask(
                    populateMolSet,
                    instance=instance,
                    eager=hasattr(settings, 'CELERY_TASK_ALWAYS_EAGER') and settings.CELERY_TASK_ALWAYS_EAGER,
                    args=(
                        instance.pk,
                        self.get_initializer_class(),
                        self.get_initializer_additional_arguments(serializer.validated_data)
                    ),
                )
                ret = serializer_class(instance).data
                ret["taskID"] = task_id
                return Response(ret, status=status.HTTP_201_CREATED)
            except Exception as exp:
                traceback.print_exc()
                if task and task.id:
                    settings.CURRENT_CELERY_APP.control.revoke(task_id=task.id, terminate=True)
                instance.delete()
                return Response({"error" : repr(exp)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
        else:
            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def update(self, request, *args, **kwargs):
        serializer_class = self.get_serializer_class()
        serializer = serializer_class(data=request.data)
        if serializer.is_valid():
            with transaction.atomic():
                instance = serializer.update(MolSet.objects.get(pk=kwargs['pk']), serializer.validated_data)

            task = None
            try:
                task, task_id = runTask(
                    updateMolSet,
                    instance=instance,
                    eager=hasattr(settings, 'CELERY_TASK_ALWAYS_EAGER') and settings.CELERY_TASK_ALWAYS_EAGER,
                    args=(
                        instance.pk,
                        self.get_updater_class(),
                        self.get_updater_additional_arguments(serializer.validated_data)
                    ),
                )
                ret = serializer_class(instance).data
                ret["taskID"] = task_id
                return Response(ret, status=status.HTTP_202_ACCEPTED)
            except Exception as exp:
                traceback.print_exc()
                if task and task.id:
                    settings.CURRENT_CELERY_APP.control.revoke(task_id=task.id, terminate=True)
                return Response({"error" : repr(exp)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

class ActivitySetViewSet(
    FilterToProjectMixIn,
    FilterToUserMixIn,
    mixins.ListModelMixin,
    mixins.RetrieveModelMixin,
    mixins.DestroyModelMixin,
    viewsets.GenericViewSet
):
    queryset = ActivitySet.objects.all()
    serializer_class = ActivitySetSerializer
    owner_relation = "project__owner"

    project_id_param = openapi.Parameter('project_id', openapi.IN_QUERY, description="Return activity sets related to just this project.", type=openapi.TYPE_NUMBER)
    @swagger_auto_schema(
        operation_description="List all activity sets. Can give a project ID to filter on."
        # , methods=['GET']
        , manual_parameters=[project_id_param]
        , responses={200: ActivitySetSerializer(many=True)}
    )
    def list(self, request, *args, **kwargs):
        return super().list(request, *args, **kwargs)

    mols_param = openapi.Parameter('mols', openapi.IN_QUERY, description="Only return activities for the given molecule IDs.", type=openapi.TYPE_ARRAY, items=openapi.Items(type=openapi.TYPE_INTEGER))
    @swagger_auto_schema(
        methods=['GET']
        , manual_parameters=[mols_param]
        , responses={200: ActivitySerializer(many=True)}
    )
    @action(detail=True, methods=['get'])
    def activities(self, request, pk=None):
        try:
            activity_set = self.get_queryset().get(pk=pk)
        except ActivitySet.DoesNotExist:
            return Response({"error" : f"No such set: {pk}"}, status=status.HTTP_404_NOT_FOUND)
        activities = Activity.objects.filter(source=activity_set)

        mols = self.request.query_params.get('mols', [])
        if mols:
            mols = mols.split(',')
            activities = activities.filter(molecule__in=mols)

        activities =  activities.order_by('id')
        paginator = ActivityPagination()
        page = paginator.paginate_queryset(activities, self.request, view=self)
        if page is not None:
            serializer = ActivitySerializer(page, many=True)
            return paginator.get_paginated_response(serializer.data)
        else:
            return Response({"error" : "You need to specify a valid page number."}, status=status.HTTP_400_BAD_REQUEST)

    @swagger_auto_schema(
        methods=['GET']
        , responses={200: ActivitySetSummarySerializer(many=False)}
    )
    @action(detail=True, methods=['get'])
    def summary(self, request, pk):
        try:
            activity_set = self.get_queryset().get(pk=pk)
        except ActivitySet.DoesNotExist:
            return Response({"error" : f"No such set: {pk}"}, status=status.HTTP_404_NOT_FOUND)

        summary = activity_set.getSummary()
        serializer = ActivitySetSummarySerializer(summary)
        data = serializer.data
        return Response(data, status=status.HTTP_200_OK)

class GeneratedSetViewSet(BaseMolSetViewSet):
    queryset = GeneratedMolSet.objects.all()
    serializer_class = GeneratedSetSerializer
    initializer_class = GeneratedSetInitializer

    def get_serializer_class(self):
        if self.action in ('create', 'update', 'partial_update'):
            return GeneratedSetInitSerializer
        else:
            return super().get_serializer_class()

    def get_initializer_additional_arguments(self, validated_data):
        return {
            "n_samples" : validated_data["nSamples"]
        }

class MolSetMoleculesView(generics.ListAPIView):
    pagination_class = MoleculePagination
    queryset = Molecule.objects.order_by('id')

    # FIXME: this action is paginated, but it needs to be indicated in the swagger docs somehow
    @swagger_auto_schema(responses={200: MoleculeSerializer(many=True)})
    def get(self, request, pk):
        try:
            molset = MolSet.objects.get(pk=pk)
            if molset.project.owner != request.user:
                raise MolSet.DoesNotExist
        except MolSet.DoesNotExist:
            return Response({"error" : f"No such set: {pk}"}, status=status.HTTP_404_NOT_FOUND)
        molset_mols = self.get_queryset().filter(providers__id = molset.id)
        page = self.paginate_queryset(molset_mols)
        if page is not None:
            serializer = MoleculeSerializer(page, many=True, context={"request": request})
            return self.get_paginated_response(serializer.data)
        else:
            return Response({"error" : "You need to specify a valid page number."}, status=status.HTTP_400_BAD_REQUEST)

class MoleculeViewSet(
                   FilterToUserMixIn,
                   mixins.RetrieveModelMixin,
                   mixins.DestroyModelMixin,
                   GenericViewSet):
    queryset = Molecule.objects.order_by('id')
    serializer_class = MoleculeSerializer
    pagination_class = MoleculePagination
    owner_relation = 'providers__project__owner'

    def get_queryset(self):
        ret = super().get_queryset().distinct()

        if self.action in ('retrieve', 'list') and 'properties' in self.request.query_params:
            for prop in self.request.query_params['properties'].split(','):
                lookup = f"rdkit_prop_{prop}"
                prop_calculator = getattr(djrdkit, prop)
                ret = ret.annotate(**{ lookup: prop_calculator('molObject')})
        return ret

    properties = openapi.Parameter('properties', openapi.IN_QUERY, description="Attach specified physchem properties to the response. You should be able to use all properties listed here: https://django-rdkit.readthedocs.io/en/latest/functions.html", type=openapi.TYPE_ARRAY, items=openapi.Items(type=openapi.TYPE_STRING))
    @swagger_auto_schema(
        manual_parameters=[properties]
    )
    def retrieve(self, request, *args, **kwargs):
        return super().retrieve(request, *args, **kwargs)

    # FIXME: this action is paginated, but it needs to be indicated in the swagger docs somehow
    molset_id_param = openapi.Parameter('activity_set', openapi.IN_QUERY, description="Return only activities that belong to a certain activity set.", type=openapi.TYPE_NUMBER)
    @swagger_auto_schema(
        methods=['GET']
        , responses={200: ActivitySerializer(many=True)}
        , manual_parameters=[molset_id_param]
    )
    @action(detail=True, methods=['get'])
    def activities(self, request, pk=None):
        try:
            molecule = self.get_queryset().get(pk=pk)
        except Molecule.DoesNotExist:
            return Response({"error" : f"No molecule found: {pk}"}, status=status.HTTP_404_NOT_FOUND)
        activities = molecule.activities.filter(source__project__owner=self.request.user).distinct().order_by('id')
        activity_set = self.request.query_params.get('activity_set', None)
        if activity_set is not None:
            activities = activities.filter(source__pk=int(activity_set))
        paginator = ActivityPagination()
        page = paginator.paginate_queryset(activities, self.request, view=self)
        if page is not None:
            serializer = ActivitySerializer(page, many=True)
            return paginator.get_paginated_response(serializer.data)
        else:
            return Response({"error" : "You need to specify a valid page number."}, status=status.HTTP_400_BAD_REQUEST)

class MolSetViewSet(
    FilterToProjectMixIn
    , FilterToUserMixIn
    , mixins.ListModelMixin
    , mixins.DestroyModelMixin
    , GenericViewSet
):
    queryset = MolSet.objects.order_by('id')
    serializer_class = GenericMolSetSerializer
    owner_relation = "project__owner"

    project_id_param = openapi.Parameter('project_id', openapi.IN_QUERY, description="Return compound sets related to just this project.", type=openapi.TYPE_NUMBER)
    @swagger_auto_schema(
        operation_description="List all compound sets. Can give a project ID to filter on."
        # , methods=['GET']
        , manual_parameters=[project_id_param]
        , responses={200: GenericMolSetSerializer(many=True)}
    )
    def list(self, request, *args, **kwargs):
        return super().list(request, *args, **kwargs)
