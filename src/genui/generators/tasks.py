"""
tasks

Created by: Martin Sicho
On: 28-01-20, 13:52
"""
from celery import shared_task

from commons.tasks import ProgressRecorder
from generators import models

from generators.core import builders


@shared_task(name="buildDrugExModel", bind=True)
def buildDrugExModel(self, model_id, builder_class, model_class):
    model_class = getattr(models, model_class)
    instance = model_class.objects.get(pk=model_id)
    builder_class = getattr(builders, builder_class)
    recorder = ProgressRecorder(self)
    if hasattr(instance, 'parent'):
        builder = builder_class(
            instance,
            instance.parent,
            progress=recorder
        )
    else:
        builder = builder_class(
            instance,
            progress=recorder
        )
    builder.build()

    return {
        "errors" : [repr(x) for x in builder.errors],
        "modelFile" : instance.modelFile.path
    }
