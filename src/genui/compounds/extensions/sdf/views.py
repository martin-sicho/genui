from rest_framework.parsers import MultiPartParser, JSONParser

from . import models
from . import serializers
from .initializer import SDFSetInitializer
from genui.compounds.views import BaseMolSetViewSet
from . import apps

class SDFSetViewSet(BaseMolSetViewSet):
    parser_classes = (MultiPartParser, JSONParser)
    queryset = models.SDFCompounds.objects.all()
    serializer_class = serializers.SDFSetSerializer
    initializer_class = SDFSetInitializer

    def get_serializer_class(self):
        if self.action in ('update', 'partial_update',):
            return serializers.SDFSetUpdateSerializer
        else:
            return super().get_serializer_class()

    def get_initializer_additional_arguments(self, validated_data):
        return {
            'parser_class' : f'{apps.SdfConfig.name}.parser.SDFParser'
        }