import traceback

from django.conf import settings
from django.db import transaction
from drf_yasg import openapi
from drf_yasg.utils import swagger_auto_schema
from rest_framework import viewsets, mixins, status
from rest_framework.response import Response

from commons.views import FilterToProjectMixIn
from qsar.tasks import buildModel
from . import models
from . import serializers


class QSARModelViewSet(FilterToProjectMixIn, viewsets.ModelViewSet):
    queryset = models.QSARModel.objects.all()
    serializer_class = serializers.QSARModelSerializer

    def get_serializer_class(self):
        if self.action == 'create':
            return serializers.QSARModelInitSerializer
        # elif self.action in ('update', 'partial_update',):
        #     return ChEMBLSetUpdateSerializer
        else:
            return super().get_serializer_class()

    project_id_param = openapi.Parameter('project_id', openapi.IN_QUERY, description="ID of a project to limit the list of results to.", type=openapi.TYPE_NUMBER)
    @swagger_auto_schema(
        operation_description="List all models. Supply a project ID to get only models specific to a particular project."
        # , methods=['GET']
        , manual_parameters=[project_id_param]
        # , responses={200: GenericMolSetSerializer(many=True)}
    )
    def list(self, request, *args, **kwargs):
        return super().list(request, *args, **kwargs)

    def create(self, request, *args, **kwargs):
        serializer = serializers.QSARModelInitSerializer(data=request.data)
        if serializer.is_valid():
            with transaction.atomic():
                instance = serializer.create(serializer.validated_data)

            task = None
            try:
                task = instance.apply_async(buildModel, args=[instance.pk, 'BasicQSARModelBuilder'])
                ret = serializers.QSARModelSerializer(instance).data
                ret["taskID"] = task.id
                return Response(ret, status=status.HTTP_201_CREATED)
            except Exception as exp:
                traceback.print_exc()
                if task and task.id:
                    settings.CURRENT_CELERY_APP.control.revoke(task_id=task.id, terminate=True)
                instance.delete()
                return Response({"error" : repr(exp)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
        print(serializer.errors)
        print(serializer.initial_data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


class DescriptorGroupsViewSet(
    mixins.ListModelMixin,
    mixins.RetrieveModelMixin,
    viewsets.GenericViewSet
):
    queryset = models.DescriptorGroup.objects.all()
    serializer_class = serializers.DescriptorGroupSerializer

