"""
views

Created by: Martin Sicho
On: 1/1/20, 4:32 PM
"""

from drf_yasg.utils import swagger_auto_schema
from rest_framework import views, status
from celery_progress.backend import Progress
from rest_framework.response import Response

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
        return Response({"error" : "Invalid serializer instance. Did you provide a valid task ID?"}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

