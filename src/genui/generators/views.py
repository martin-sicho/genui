from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema
from rest_framework import mixins
from rest_framework.viewsets import GenericViewSet

from genui.accounts.serializers import FilterToUserMixIn
from genui.projects.serializers import FilterToProjectMixIn
from genui.modelling.core.bases import Algorithm
from genui.modelling.models import AlgorithmMode
from genui.modelling.views import ModelViewSet, AlgorithmViewSet, MetricsViewSet
from . import models
from . import serializers
from .core import builders
from .core import algorithms
from .tasks import buildDrugExModel

class GeneratorViewSet(
        FilterToProjectMixIn
        , FilterToUserMixIn
        , mixins.ListModelMixin
        , mixins.DestroyModelMixin
        , GenericViewSet
    ):
    queryset = models.Generator.objects.all()
    serializer_class = serializers.GeneratorSerializer
    owner_relation = 'project__owner'

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
    build_task = buildDrugExModel

    def get_builder_kwargs(self):
        return {"model_class" : models.DrugExNet.__name__}

class DrugExAgentViewSet(ModelViewSet):
    queryset = models.DrugExAgent.objects.all()
    serializer_class = serializers.DrugExAgentSerializer
    init_serializer_class = serializers.DrugExAgentInitSerializer
    builder_class = builders.DrugExAgentBuilder
    build_task = buildDrugExModel

    def get_builder_kwargs(self):
        return {"model_class" : models.DrugExAgent.__name__}

class GeneratorAlgorithmViewSet(AlgorithmViewSet):

    def get_queryset(self):
        current = super().get_queryset()
        return current.filter(validModes__name__in=(algorithms.DrugExNetwork.GENERATOR,)).distinct('id')

class GeneratorMetricsViewSet(MetricsViewSet):

    def get_queryset(self):
        ret = super().get_queryset()
        modes = AlgorithmMode.objects.filter(name__in=(Algorithm.GENERATOR,))
        return ret.filter(validModes__in=modes)