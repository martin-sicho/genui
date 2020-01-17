"""
views

Created by: Martin Sicho
On: 1/1/20, 4:32 PM
"""

from drf_yasg.utils import swagger_auto_schema
from rest_framework import views, status
from celery_progress.backend import Progress
from rest_framework.response import Response
from rest_framework.schemas.openapi import AutoSchema

from commons.serializers import TasksSerializerFactory
from compounds.models import MolSet

from .serializers import TaskProgressSerializer


class TaskProgressView(views.APIView):

    @swagger_auto_schema(responses={200: TaskProgressSerializer()})
    def get(self, request, task_id):
        progress = Progress(task_id)
        try:
            info = progress.get_info()
        except Exception as exp:
            return Response({"error" : "Failed to retrieve progress info. The following exception occurred: " + repr(exp)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
        serializer = TaskProgressSerializer(data=info)
        if serializer.is_valid():
            return Response(info, status=status.HTTP_200_OK)
        ret = serializer.errors
        ret.update({"error" : "Invalid serializer instance. Did you provide a valid task ID?"})
        return Response(ret, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


class ModelTasksView(views.APIView):
    class Schema(TasksSerializerFactory.AutoSchemaMixIn, AutoSchema):
        pass

    started_only = False
    model_class = None # needs to implement getTasksAsDict
    schema = Schema()

    @swagger_auto_schema(responses={200: TasksSerializerFactory.get(["someTaskName"])})
    def get(self, request, pk):
        try:
            molset = self.model_class.objects.get(pk=pk)
        except MolSet.DoesNotExist:
            return Response({"error" : f"No such set. Unknown ID: {pk}"}, status=status.HTTP_400_BAD_REQUEST)
        data = molset.getTasksAsDict(self.started_only)
        ser = TasksSerializerFactory.get(data.keys())
        serializer = ser(
            data=data
        )
        if serializer.is_valid():
            return Response(serializer.data, status=status.HTTP_200_OK)
        return Response(serializer.errors, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

class FilterToProjectMixIn:

    def get_queryset(self):
        queryset = super().get_queryset()
        project = self.request.query_params.get('project_id', None)
        if project is not None:
            queryset = queryset.filter(project__pk=int(project))
        return queryset