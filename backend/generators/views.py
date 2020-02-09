from django.shortcuts import render
from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema
from rest_framework import mixins
from rest_framework.viewsets import GenericViewSet

from commons.views import FilterToProjectMixIn
from modelling.views import ModelViewSet, AlgorithmViewSet
from . import models
from . import serializers
from .core import builders
from .core import algorithms
from .tasks import buildDrugExNet, buildDrugExAgent

class GeneratorViewSet(
        FilterToProjectMixIn
        , mixins.ListModelMixin
        , mixins.DestroyModelMixin
        , GenericViewSet
    ):
    queryset = models.Generator.objects.all()
    serializer_class = serializers.GeneratorSerializer

    project_id_param = openapi.Parameter('project_id', openapi.IN_QUERY, description="Return generators related to just the project with this ID.", type=openapi.TYPE_NUMBER)
    @swagger_auto_schema(
        operation_description="List all available generators. Can give a project ID to filter on."
        # , methods=['GET']
        , manual_parameters=[project_id_param]
        , responses={200: serializers.GeneratorSerializer(many=True)}
    )
    def list(self, request, *args, **kwargs):
        return super().list(request, *args, **kwargs)

class DrugExNetViewSet(ModelViewSet):
    queryset = models.DrugExNet.objects.all()
    serializer_class = serializers.DrugExNetSerializer
    init_serializer_class = serializers.DrugExNetInitSerializer
    builder_class = builders.DrugExNetBuilder
    build_task = buildDrugExNet

class DrugExAgentViewSet(ModelViewSet):
    queryset = models.DrugExAgent.objects.all()
    serializer_class = serializers.DrugExAgentSerializer
    init_serializer_class = serializers.DrugExAgentInitSerializer
    builder_class = builders.DrugExAgentBuilder
    build_task = buildDrugExAgent

class GeneratorAlgorithmViewSet(AlgorithmViewSet):

    def get_queryset(self):
        current = super().get_queryset()
        return current.filter(validModes__name__in=(algorithms.DrugExNetwork.GENERATOR,)).distinct('id')