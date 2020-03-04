"""
tasks

Created by: Martin Sicho
On: 25-02-20, 16:37
"""
from celery import shared_task

from commons.tasks import ProgressRecorder
from maps.models import Map
from .core import builders


@shared_task(name="CreateMap", bind=True)
def createMap(self, model_id, builder_class):
    instance = Map.objects.get(pk=model_id)
    builder_class = getattr(builders, builder_class)
    recorder = ProgressRecorder(self)
    builder = builder_class(
        instance,
        recorder
    )
    builder.build()

    return {
        "errors" : [repr(x) for x in builder.errors],
        "mapFile" : instance.modelFile.path
    }

