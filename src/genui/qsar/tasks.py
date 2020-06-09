"""
tasks

Created by: Martin Sicho
On: 29-11-19, 13:44
"""

from genui.utils.extensions.tasks.progress import ProgressRecorder
from celery import shared_task

from .models import QSARModel, ModelActivitySet
from genui.utils.inspection import getObjectAndModuleFromFullName


@shared_task(name="BuildQSARModel", bind=True)
def buildQSARModel(self, model_id, builder_class):
    instance = QSARModel.objects.get(pk=model_id)
    builder_class = getObjectAndModuleFromFullName(builder_class)[0]
    recorder = ProgressRecorder(self)
    builder = builder_class(
        instance,
        recorder
    )
    builder.build()

    return {
        "errors" : [repr(x) for x in builder.errors],
        "modelName" : instance.name,
        "modelID" : instance.id,
    }

@shared_task(name="PredictWithQSARModel", bind=True)
def predictWithQSARModel(self, predictions_id, builder_class):
    instance = ModelActivitySet.objects.get(pk=predictions_id)
    model = QSARModel.objects.get(pk=instance.model.id)
    builder_class = getObjectAndModuleFromFullName(builder_class)[0]
    recorder = ProgressRecorder(self)
    builder = builder_class(
        model,
        recorder
    )
    builder.populateActivitySet(instance)

    return {
        "errors" : [repr(x) for x in builder.errors],
        "modelName" : model.name,
        "modelID" : model.id,
        "activitySetName" : instance.name,
        "activitySetID" : instance.id,
    }