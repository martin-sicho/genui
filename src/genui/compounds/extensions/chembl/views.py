from rest_framework import mixins, viewsets

from genui.compounds.views import BaseMolSetViewSet
from . import models
from . import serializers
from .initializers import ChEMBLSetInitializer

class ChEMBLSetViewSet(BaseMolSetViewSet):
    queryset = models.ChEMBLCompounds.objects.all()
    serializer_class = serializers.ChEMBLSetSerializer
    initializer_class = ChEMBLSetInitializer

    def get_serializer_class(self):
        if self.action == 'create':
            return serializers.ChEMBLSetInitSerializer
        elif self.action in ('update', 'partial_update',):
            return serializers.ChEMBLSetUpdateSerializer
        else:
            return super().get_serializer_class()

    def get_initializer_additional_arguments(self, validated_data):
        return {
                    "targets" : list(set(validated_data["targets"])),
                    "max_per_target" : validated_data["maxPerTarget"] if "maxPerTarget" in validated_data else None
        }

class ChEMBLAssayViewSet(
    mixins.RetrieveModelMixin,
    viewsets.GenericViewSet
):
    queryset = models.ChEMBLAssay.objects.all()
    serializer_class = serializers.ChEMBLAssaySerializer

class ChEMBLTargetViewSet(
    mixins.RetrieveModelMixin,
    viewsets.GenericViewSet
):
    queryset = models.ChEMBLTarget.objects.all()
    serializer_class = serializers.ChEMBLTargetSerializer
