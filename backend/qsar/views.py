import traceback

from django.conf import settings
from rest_framework import viewsets, mixins, status
from rest_framework.response import Response

from modelling.views import ModelViewSet, AlgorithmViewSet
from qsar.core.builders import BasicQSARModelBuilder
from qsar.core.bases import Algorithm
from . import models
from . import serializers
from .tasks import buildModel, predictWithModel


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

class ModelPredictionsViewSet(viewsets.ModelViewSet):
    queryset = models.ModelActivitySet.objects.all()
    serializer_class = serializers.ModelActivitySetSerializer

    def create(self, request, *args, **kwargs):
        created = super().create(request, *args, **kwargs)


        task = None
        try:
            instance = models.QSARModel.objects.get(pk=created.data['model'])
            task = instance.apply_async(predictWithModel, args=[instance.pk, instance.builder.name])
            created.data["taskID"] = task.id
            return created
        except Exception as exp:
            traceback.print_exc()
            if task and task.id:
                settings.CURRENT_CELERY_APP.control.revoke(task_id=task.id, terminate=True)

            models.ModelActivitySet.objects.get(pk=created.data['id']).delete()
            return Response({"error" : repr(exp)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)