from django.shortcuts import render

from modelling.views import ModelViewSet, AlgorithmViewSet
from . import models
from . import serializers
from .core import builders
from .core import algorithms
from .tasks import buildDrugExNet

class DrugExNetViewSet(ModelViewSet):
    queryset = models.DrugExNet.objects.all()
    serializer_class = serializers.DrugExNetSerializer
    init_serializer_class = serializers.DrugExNetInitSerializer
    builder_class = builders.DrugExBuilder
    build_task = buildDrugExNet

class DrugExNetAlgorithmViewSet(AlgorithmViewSet):

    def get_queryset(self):
        current = super().get_queryset()
        return current.filter(validModes__name__in=(algorithms.DrugExNetwork.GENERATOR,)).distinct('id')