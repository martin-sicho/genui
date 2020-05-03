import traceback

from django.conf import settings
from drf_yasg.utils import swagger_auto_schema
from rest_framework import viewsets, mixins, status
from rest_framework.decorators import action
from rest_framework.response import Response

from genui.modelling.models import AlgorithmMode
from genui.modelling.views import ModelViewSet, AlgorithmViewSet, MetricsViewSet, PredictMixIn
from genui.qsar.genuimodels.builders import BasicQSARModelBuilder
from genui.qsar.genuimodels.bases import Algorithm
from . import models
from . import serializers
from .tasks import buildQSARModel, predictWithQSARModel
from genui.utils.extensions.tasks.utils import runTask


class QSARModelViewSet(PredictMixIn, ModelViewSet):
    queryset = models.QSARModel.objects.all()
    serializer_class = serializers.QSARModelSerializer
    init_serializer_class = serializers.QSARModelInitSerializer
    builder_class = BasicQSARModelBuilder
    build_task = buildQSARModel
    predict_task = predictWithQSARModel

    @swagger_auto_schema(
        methods=['GET']
        , responses={
            200: serializers.ModelActivitySetSerializer(many=True),
        }
    )
    @swagger_auto_schema(
        methods=['POST']
        , responses={
            201: serializers.ModelActivitySetSerializer(many=False),
        }
        , request_body=serializers.ModelActivitySetSerializer(many=False)
    )
    @swagger_auto_schema(
        methods=['DELETE']
        , responses={
            204: "",
        }
    )
    @action(detail=True, methods=['get', 'post', 'delete'])
    def predictions(self, request, pk=None):
        # FIXME: some of this should be moved to the PredictMixIn for reuse
        try:
            instance = self.get_queryset().get(pk=pk)
        except models.QSARModel.DoesNotExist:
            return Response({"error" : f"Model not found: {pk}"}, status=status.HTTP_404_NOT_FOUND)

        if request.method == 'GET':
            predictions = instance.predictions.all()
            serializer = serializers.ModelActivitySetSerializer(predictions, many=True)
            return Response(serializer.data)

        elif request.method == 'POST':
            request.data['project'] = instance.project.id
            request.data['model'] = instance.id
            serializer = serializers.ModelActivitySetSerializer(data=request.data, many=False)
            if serializer.is_valid():
                created = serializer.create(serializer.validated_data)

                task = None
                try:
                    task, task_id = runTask(
                        self.get_predict_task(),
                        instance=instance,
                        eager=hasattr(settings, 'CELERY_TASK_ALWAYS_EAGER') and settings.CELERY_TASK_ALWAYS_EAGER,
                        args=(
                            created.pk,
                            self.get_builder_class()
                        )
                    )
                    ret = serializers.ModelActivitySetSerializer(created, many=False)
                    ret.data["taskID"] = task_id
                    return Response(ret.data, status=status.HTTP_201_CREATED)
                except Exception as exp:
                    traceback.print_exc()
                    if task and task.id:
                        settings.CURRENT_CELERY_APP.control.revoke(task_id=task.id, terminate=True)
                    created.delete()
                    return Response({"error" : repr(exp)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

        elif request.method == 'DELETE':
            instance.predictions.all().delete()
            return Response(status=status.HTTP_204_NO_CONTENT)


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

class QSARMetricsViewSet(MetricsViewSet):

    def get_queryset(self):
        ret = super().get_queryset()
        modes = AlgorithmMode.objects.filter(name__in=(Algorithm.CLASSIFICATION, Algorithm.REGRESSION))
        return ret.filter(validModes__in=modes)