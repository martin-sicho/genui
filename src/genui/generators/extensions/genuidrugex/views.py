import traceback

from django.conf import settings
from rest_framework import viewsets, status
from rest_framework.decorators import action
from rest_framework.response import Response

from genui import celery_app
from genui.accounts.serializers import FilterToUserMixIn
from genui.projects.serializers import FilterToProjectMixIn
from genui.utils.extensions.tasks.utils import runTask
from . import models
from . import serializers
from .genuimodels import builders
from .tasks import buildDrugExModel, calculateEnvironment
from genui.models.views import ModelViewSet


class DrugExNetViewSet(ModelViewSet):
    queryset = models.DrugExNet.objects.order_by('-created')
    serializer_class = serializers.DrugExNetSerializer
    init_serializer_class = serializers.DrugExNetInitSerializer
    builder_class = builders.DrugExNetBuilder
    build_task = buildDrugExModel

    def get_builder_kwargs(self):
        return {"model_class" : models.DrugExNet.__name__}


class DrugExAgentViewSet(ModelViewSet):
    queryset = models.DrugExAgent.objects.order_by('-created')
    serializer_class = serializers.DrugExAgentSerializer
    init_serializer_class = serializers.DrugExAgentInitSerializer
    builder_class = builders.DrugExAgentBuilder
    build_task = buildDrugExModel

    def get_builder_kwargs(self):
        return {"model_class" : models.DrugExAgent.__name__}

class EnvironmentViewSet(FilterToProjectMixIn, FilterToUserMixIn, viewsets.ModelViewSet):
    queryset = models.DrugExEnvironment.objects.order_by('-created')
    serializer_class = serializers.DrugExEnvironmentSerializer
    calculation_serializer_class = serializers.DrugExEnvironmentCalculationSerializer
    owner_relation = "project__owner"

    @action(detail=True, methods=['post'])
    def calculate(self, request, pk=None):
        environment = self.queryset.get(pk=pk)
        serializer = self.get_serializer(data=request.data)
        if serializer.is_valid():
            data = serializer.validated_data
            task = None
            try:
                task, task_id = runTask(
                    calculateEnvironment,
                    eager=hasattr(settings, 'CELERY_TASK_ALWAYS_EAGER') and settings.CELERY_TASK_ALWAYS_EAGER,
                    args=(
                        environment.pk,
                        data["molsets"],
                        data["useModifiers"]
                    ),
                )
                data["taskID"] = task_id
                return Response(data, status=status.HTTP_201_CREATED)
            except Exception as exp:
                traceback.print_exc()
                if task and task.id:
                    celery_app.control.revoke(task_id=task.id, terminate=True)
                return Response({"error" : repr(exp)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
        else:
            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def get_serializer_class(self):
        if self.action == 'calculate':
            return self.calculation_serializer_class
        else:
            return self.serializer_class

class ScoringMethodViewSet(FilterToProjectMixIn, FilterToUserMixIn, viewsets.ModelViewSet):
    queryset = models.ScoringMethod.objects.order_by('-created')
    serializer_class = serializers.ScoringFunctionSerializer
    owner_relation = "project__owner"
    # http_method_names = ['get']

class QSARScorerViewSet(ScoringMethodViewSet):
    queryset = models.GenUIModelScorer.objects.order_by('-created')
    serializer_class = serializers.QSARScorerSerializer

class PropertyScorerViewSet(ScoringMethodViewSet):
    queryset = models.PropertyScorer.objects.order_by('-created')
    serializer_class = serializers.PropertyScorerSerializer

class ModifierViewSet(FilterToProjectMixIn, FilterToUserMixIn, viewsets.ModelViewSet):
    queryset = models.ScoreModifier.objects.order_by('-created')
    serializer_class = serializers.ModifierSerializer
    test_serializer_class = serializers.ModifierTestSerializer
    owner_relation = "project__owner"

    @action(detail=False, methods=['post'])
    def test(self, request):
        modifier = self.queryset.model
        serializer = self.get_serializer(data=request.data)
        if serializer.is_valid():
            data = serializer.validated_data
            results = modifier.test(data['inputs'], **(data['params']))
            data["results"] = results
            return Response(data)
        else:
            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def get_serializer_class(self):
        if self.action == 'test':
            return self.test_serializer_class
        else:
            return self.serializer_class

class ClippedViewSet(ModifierViewSet):
    queryset = models.ClippedScore.objects.order_by('-created')
    serializer_class = serializers.ClippedSerializer

class HumpViewSet(ModifierViewSet):
    queryset = models.SmoothHump.objects.order_by('-created')
    serializer_class = serializers.SmoothHumpSerializer

class ScorerViewSet(FilterToProjectMixIn, FilterToUserMixIn, viewsets.ModelViewSet):
    queryset = models.DrugExScorer.objects.order_by('-created')
    serializer_class = serializers.DrugExScorerSerializer
    owner_relation = "project__owner"

class GeneratorViewSet(FilterToProjectMixIn, FilterToUserMixIn, viewsets.ModelViewSet):
    queryset = models.DrugEx.objects.order_by('-created')
    serializer_class = serializers.DrugExGeneratorSerializer
    owner_relation = "project__owner"