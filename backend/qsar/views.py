from rest_framework import viewsets, mixins

from modelling.views import ModelViewSet, AlgorithmViewSet
from qsar.core.builders import BasicQSARModelBuilder
from qsar.core.bases import Algorithm
from . import models
from . import serializers
from .tasks import buildModel

class QSARModelViewSet(ModelViewSet):
    queryset = models.QSARModel.objects.all()
    serializer_class = serializers.QSARModelSerializer
    init_serializer_class = serializers.QSARModelInitSerializer
    builder_class = BasicQSARModelBuilder
    build_task = buildModel


class DescriptorGroupsViewSet(
    mixins.ListModelMixin,
    mixins.RetrieveModelMixin,
    viewsets.GenericViewSet
):
    queryset = models.DescriptorGroup.objects.all()
    serializer_class = serializers.DescriptorGroupSerializer

class QSARAlgorithmViewSet(AlgorithmViewSet):

    def get_queryset(self):
        current = super().get_queryset()
        return current.filter(validModes__name__in=(Algorithm.CLASSIFICATION, Algorithm.REGRESSION)).distinct('id')

