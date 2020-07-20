from rest_framework.parsers import MultiPartParser

from genui.compounds.views import BaseMolSetViewSet
from . import models
from . import serializers
from . import apps
from .initializer import CSVSetInitializer

class CSVSetViewSet(BaseMolSetViewSet):
    parser_classes = (MultiPartParser,)
    queryset = models.CSVCompounds.objects.all()
    serializer_class = serializers.CSVSetSerializer
    initializer_class = CSVSetInitializer

    def get_initializer_additional_arguments(self, validated_data):
        return {
            'parser_class' : f'{apps.CsvimportsConfig.name}.parser.CSVParser'
        }
