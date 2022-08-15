"""
tasks

Created by: Martin Sicho
On: 28-01-20, 13:52
"""
from celery import shared_task

from genui.utils.extensions.tasks.progress import ProgressRecorder
from genui.utils.inspection import getObjectAndModuleFromFullName

from . import models
from .torchutils import cleanup

@shared_task(name="BuildDrugExModel", bind=True, queue='gpu')
def buildDrugExModel(self, model_id, builder_class, model_class):
    # get the builder
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

    # build the model
    try:
        builder.build()
    except Exception as exp:
        raise exp

    cleanup()
    return {
        "errors" : [repr(x) for x in builder.errors],
        "DrExModelName" : instance.name,
        "DrExModelID" : instance.id,
    }
