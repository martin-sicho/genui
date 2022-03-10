from celery.result import AsyncResult
from celery_progress.backend import Progress

# Create your views here.
from drf_yasg.utils import swagger_auto_schema
from rest_framework import permissions, views, status
from rest_framework.decorators import permission_classes
from rest_framework.response import Response
from rest_framework.schemas.openapi import AutoSchema

from .serializers import TaskProgressSerializer, TasksSerializerFactory


@permission_classes((permissions.IsAuthenticated,))
class TaskProgressView(views.APIView):

    @swagger_auto_schema(responses={200: TaskProgressSerializer()})
    def get(self, request, task_id):
        # TODO: check if the task belongs to the logged in user
        result = AsyncResult(task_id)
        progress = Progress(result)
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
    permission_classes = (permissions.IsAuthenticated,)

    @swagger_auto_schema(responses={200: TasksSerializerFactory.get(["someTaskName"])})
    def get(self, request, pk):
        try:
            model = self.model_class.objects.get(pk=pk)
            if hasattr(model, 'project') and (model.project.owner != request.user):
                raise self.model_class.DoesNotExist
        except self.model_class.DoesNotExist:
            return Response({"error" : f"No model found: {pk}"}, status=status.HTTP_404_NOT_FOUND)
        data = model.getTasksAsDict(self.started_only)
        ser = TasksSerializerFactory.get(data.keys())
        serializer = ser(
            data=data
        )
        if serializer.is_valid():
            return Response(serializer.data, status=status.HTTP_200_OK)
        return Response(serializer.errors, status=status.HTTP_500_INTERNAL_SERVER_ERROR)