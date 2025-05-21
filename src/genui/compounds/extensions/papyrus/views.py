from rest_framework import mixins, viewsets

from genui.compounds.views import BaseMolSetViewSet
from . import models
from . import serializers
from .initializers import PapyrusSetInitializer

class PapyrusSetViewSet(BaseMolSetViewSet):
    queryset = models.PapyrusCompounds.objects.all()
    serializer_class = serializers.PapyrusSetSerializer
    initializer_class = PapyrusSetInitializer

    def get_serializer_class(self):
        if self.action == 'create':
            return serializers.PapyrusSetInitSerializer
        elif self.action in ('update', 'partial_update',):
            return serializers.PapyrusSetUpdateSerializer
        else:
            return super().get_serializer_class()

    def get_initializer_additional_arguments(self, validated_data):
        return {
                    "targets" : list(set(validated_data["targets"])),
                    "max_per_target" : validated_data["maxPerTarget"] if "maxPerTarget" in validated_data else None
        }

class PapyrusAssayViewSet(
    mixins.RetrieveModelMixin,
    viewsets.GenericViewSet
):
    queryset = models.PapyrusAssay.objects.all()
    serializer_class = serializers.PapyrusAssaySerializer

class PapyrusTargetViewSet(
    mixins.RetrieveModelMixin,
    viewsets.GenericViewSet
):
    queryset = models.PapyrusTarget.objects.all()
    serializer_class = serializers.PapyrusTargetSerializer
