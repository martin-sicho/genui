from django.shortcuts import render

# Create your views here.
from genui.compounds.extensions.generated.initializers import GeneratedSetInitializer
from genui.compounds.views import BaseMolSetViewSet
from genui.compounds.extensions.generated.models import GeneratedMolSet
from genui.compounds.extensions.generated.serializers import GeneratedSetSerializer, GeneratedSetInitSerializer, GeneratedSetUpdateSerializer


class GeneratedSetViewSet(BaseMolSetViewSet):
    queryset = GeneratedMolSet.objects.all()
    serializer_class = GeneratedSetSerializer
    initializer_class = GeneratedSetInitializer
    gpu_support = True

    def get_serializer_class(self):
        if self.action in ('create',):
            return GeneratedSetInitSerializer
        elif self.action in ('update', 'partial_update'):
            return GeneratedSetUpdateSerializer
        else:
            return super().get_serializer_class()

    def get_initializer_additional_arguments(self, validated_data):
        return {
            "n_samples" : validated_data["nSamples"]
        }