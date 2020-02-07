import traceback

from django.db import transaction
from django.conf import settings
from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema
from rest_framework import viewsets, pagination, mixins, status, generics
from rest_framework.response import Response
from rest_framework.viewsets import GenericViewSet
from rest_framework.schemas.openapi import AutoSchema

from commons.views import FilterToProjectMixIn
from compounds.initializers.generated import GeneratedSetInitializer
import generators.models
import generators.serializers
from .initializers.chembl import ChEMBLSetInitializer
from .serializers import ChEMBLSetSerializer, MoleculeSerializer, MolSetSerializer, ChEMBLSetInitSerializer, \
    GenericMolSetSerializer, ChEMBLSetUpdateSerializer
from .models import ChEMBLCompounds, Molecule, MolSet
from .tasks import populateMolSet, updateMolSet


class MoleculePagination(pagination.PageNumberPagination):
    page_size = 10

class BaseMolSetViewSet(FilterToProjectMixIn, viewsets.ModelViewSet):
    class Schema(MolSetSerializer.AutoSchemaMixIn, AutoSchema):
        pass
    schema = Schema()
    initializer_class = None
    updater_class = None

    def get_initializer_class(self):
        if not self.initializer_class:
            raise Exception("Initializer class needs to be set.")
        return self.initializer_class

    def get_updater_class(self):
        if not self.updater_class:
            return self.get_initializer_class()
        return self.updater_class

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
                arguments = self.get_initializer_additional_arguments(serializer.validated_data)
                task = instance.apply_async(populateMolSet, args=[instance.pk, self.get_initializer_class().__name__, arguments])
                ret = serializer_class(instance).data
                ret["taskID"] = task.id
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
                instance = serializer.update(ChEMBLCompounds.objects.get(pk=kwargs['pk']), serializer.validated_data)

            task = None
            try:
                arguments = self.get_updater_additional_arguments(serializer.validated_data)
                task = instance.apply_async(updateMolSet, args=[instance.pk, self.get_updater_class().__name__, arguments])
                ret = serializer_class(instance).data
                ret["taskID"] = task.id
                return Response(ret, status=status.HTTP_202_ACCEPTED)
            except Exception as exp:
                traceback.print_exc()
                if task and task.id:
                    settings.CURRENT_CELERY_APP.control.revoke(task_id=task.id, terminate=True)
                return Response({"error" : repr(exp)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

class ChEMBLSetViewSet(BaseMolSetViewSet):
    queryset = ChEMBLCompounds.objects.all()
    serializer_class = ChEMBLSetSerializer
    initializer_class = ChEMBLSetInitializer

    def get_serializer_class(self):
        if self.action == 'create':
            return ChEMBLSetInitSerializer
        elif self.action in ('update', 'partial_update',):
            return ChEMBLSetUpdateSerializer
        else:
            return super().get_serializer_class()

    def get_initializer_additional_arguments(self, validated_data):
        return {
                    "targets" : list(set(validated_data["targets"])),
                    "max_per_target" : validated_data["maxPerTarget"] if "maxPerTarget" in validated_data else None
        }

class GeneratedSetViewSet(BaseMolSetViewSet):
    queryset = generators.models.GeneratedMolSet.objects.all()
    serializer_class = generators.serializers.GeneratedSetSerializer
    initializer_class = GeneratedSetInitializer

class MolSetMoleculesView(generics.ListAPIView):
    pagination_class = MoleculePagination
    queryset = Molecule.objects.order_by('id')

    @swagger_auto_schema(responses={200: MoleculeSerializer(many=True)})
    def get(self, request, pk):
        try:
            molset = MolSet.objects.get(pk=pk)
        except MolSet.DoesNotExist:
            return Response({"error" : f"No such set. Unknown ID: {pk}"}, status=status.HTTP_400_BAD_REQUEST)
        molset_mols = self.get_queryset().filter(providers__id = molset.id)
        page = self.paginate_queryset(molset_mols)
        if page is not None:
            serializer = MoleculeSerializer(page, many=True)
            return self.get_paginated_response(serializer.data)
        else:
            return Response({"error" : "You need to specify a valid page number."}, status=status.HTTP_400_BAD_REQUEST)

class MoleculeViewSet(
                   mixins.RetrieveModelMixin,
                   mixins.DestroyModelMixin,
                   GenericViewSet):
    queryset = Molecule.objects.order_by('id')
    serializer_class = MoleculeSerializer
    pagination_class = MoleculePagination

class MolSetViewSet(
    FilterToProjectMixIn
    , mixins.ListModelMixin
    , mixins.DestroyModelMixin
    , GenericViewSet
):
    queryset = MolSet.objects.order_by('id')
    serializer_class = GenericMolSetSerializer

    project_id_param = openapi.Parameter('project_id', openapi.IN_QUERY, description="Return compound sets related to just this project.", type=openapi.TYPE_NUMBER)
    @swagger_auto_schema(
        operation_description="List all compound sets. Can give a project ID to filter on."
        # , methods=['GET']
        , manual_parameters=[project_id_param]
        , responses={200: GenericMolSetSerializer(many=True)}
    )
    def list(self, request, *args, **kwargs):
        return super().list(request, *args, **kwargs)
