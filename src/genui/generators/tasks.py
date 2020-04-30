"""
tasks

Created by: Martin Sicho
On: 28-01-20, 13:52
"""
from celery import shared_task

from genui.extensions.tasks.progress import ProgressRecorder
from genui.commons.inspection import getObjectAndModuleFromFullName

from . import models

@shared_task(name="BuildDrugExModel", bind=True)
def buildDrugExModel(self, model_id, builder_class, model_class):
    model_class = getattr(models, model_class)
    instance = model_class.objects.get(pk=model_id)
    builder_class = getObjectAndModuleFromFullName(builder_class)[0]
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
