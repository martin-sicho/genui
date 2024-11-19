import logging
import traceback
from celery import shared_task
from genui.utils.extensions.tasks.progress import ProgressRecorder
from .models import Map
from genui.utils.inspection import getObjectAndModuleFromFullName

logger = logging.getLogger(__name__)

@shared_task(name="CreateMap", bind=True)
def createMap(self, model_id, builder_class):
    try:
        instance = Map.objects.get(pk=model_id)
        builder_class = getObjectAndModuleFromFullName(builder_class)[0]
        recorder = ProgressRecorder(self)
        builder = builder_class(instance, recorder)
        builder.build()
        return {
            "errors": [repr(x) for x in builder.errors],
            "mapName": instance.name,
            "mapID": instance.id,
        }
    except Exception as e:
        error_message = f"Error creating map {model_id}: {str(e)}"
        error_traceback = traceback.format_exc()
        logger.error(error_message)
        logger.error(error_traceback)
        raise Exception(error_message)