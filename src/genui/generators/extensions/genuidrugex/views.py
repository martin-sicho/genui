from . import models
from . import serializers
from .core import builders
from .tasks import buildDrugExModel
from genui.modelling.views import ModelViewSet


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