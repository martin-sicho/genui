from django.db import transaction
from django.conf import settings
from rest_framework import viewsets, pagination, mixins, status, views
from rest_framework.response import Response
from rest_framework.viewsets import GenericViewSet
from rest_framework.schemas.openapi import AutoSchema

from .serializers import ChEMBLSetSerializer, MoleculeSerializer, MolSetSerializer
from .models import ChEMBLCompounds, Molecule, MolSet
from .tasks import populateMolSet
from commons.serializers import TasksSerializerFactory

class ChEMBLSetViewSet(viewsets.ModelViewSet):
    class Schema(MolSetSerializer.AutoSchemaMixIn, AutoSchema):
        pass

    queryset = ChEMBLCompounds.objects.all()
    serializer_class = ChEMBLSetSerializer
    schema = Schema()

    def create(self, request, *args, **kwargs):
        serializer = ChEMBLSetSerializer(data=request.data)
        if serializer.is_valid():
            with transaction.atomic():
                instance = serializer.create(serializer.validated_data)
            task = None
            try:
                task = instance.apply_async(populateMolSet, args=[instance.pk, 'ChEMBLSetInitializer'])
                data = ChEMBLSetSerializer(instance).data
                data['task'] = task.id
                return Response(data, status=status.HTTP_201_CREATED)
            except Exception as exp:
                instance.delete()
                if task and task.id:
                    settings.CURRENT_CELERY_INSTANCE.control.revoke(task_id=task.id, terminate=True)
                return Response({"error" : repr(exp)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

class MolSetTasksView(views.APIView):
    class Schema(TasksSerializerFactory.AutoSchemaMixIn, AutoSchema):
        pass

    started_only = False
    schema = Schema()

    def get(self, request, pk):
        try:
            molset = MolSet.objects.get(pk=pk)
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


class MoleculePagination(pagination.PageNumberPagination):
    page_size = 10

class MoleculeViewSet(
                   mixins.RetrieveModelMixin,
                   mixins.DestroyModelMixin,
                   GenericViewSet):
    queryset = Molecule.objects.order_by('id')
    serializer_class = MoleculeSerializer
    pagination_class = MoleculePagination
