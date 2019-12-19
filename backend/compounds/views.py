from django.db import transaction
from django.conf import settings
from rest_framework import viewsets, pagination, mixins, status
from rest_framework.response import Response
from rest_framework.viewsets import GenericViewSet

from .serializers import ChEMBLSetSerializer, MoleculeSerializer
from .models import ChEMBLCompounds, Molecule
from .tasks import populateMolSet

class ChEMBLSetViewSet(viewsets.ModelViewSet):
    queryset = ChEMBLCompounds.objects.all()
    serializer_class = ChEMBLSetSerializer

    def create(self, request, *args, **kwargs):
        serializer = ChEMBLSetSerializer(data=request.data)
        if serializer.is_valid():
            with transaction.atomic():
                instance = serializer.create(serializer.validated_data)
            task = None
            try:
                task = instance.apply_async(populateMolSet, args=[instance.pk, 'ChEMBLSetInitializer'])
                data = ChEMBLSetSerializer(instance).data
                data['populateTask'] = task.id
                return Response(data, status=status.HTTP_201_CREATED)
            except Exception as exp:
                instance.delete()
                if task and task.id:
                    settings.CURRENT_CELERY_INSTANCE.control.revoke(task_id=task.id, terminate=True)
                return Response({"exception" : repr(exp)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

class MoleculePagination(pagination.PageNumberPagination):
    page_size = 10

class MoleculeViewSet(
                   mixins.RetrieveModelMixin,
                   mixins.DestroyModelMixin,
                   GenericViewSet):
    queryset = Molecule.objects.order_by('id')
    serializer_class = MoleculeSerializer
    pagination_class = MoleculePagination
