from rest_framework.parsers import MultiPartParser

from . import models
from . import serializers
from .initializer import SDFSetInitializer
from genui.compounds.views import BaseMolSetViewSet
from . import apps

class SDFSetViewSet(BaseMolSetViewSet):
    parser_classes = (MultiPartParser,)
    queryset = models.SDFCompounds.objects.all()
    serializer_class = serializers.SDFSetSerializer
    initializer_class = SDFSetInitializer

    def get_initializer_additional_arguments(self, validated_data):
        return {
            'parser_class' : f'{apps.SdfConfig.name}.parser.SDFParser'
        }