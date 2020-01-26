"""
tasks

Created by: Martin Sicho
On: 29-11-19, 13:44
"""

from commons.tasks import ProgressRecorder
from celery import shared_task

from .models import QSARModel
from .core import builders


@shared_task(name="BuildModel", bind=True)
def buildModel(self, model_id, builder_class):
    instance = QSARModel.objects.get(pk=model_id)
    builder_class = getattr(builders, builder_class)
    recorder = ProgressRecorder(self)
    builder = builder_class(
        instance,
        recorder
    )
    builder.build()

    return {
        "errors" : [repr(x) for x in builder.errors],
        "modelFile" : instance.modelFile.path
    }