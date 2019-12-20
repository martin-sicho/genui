"""
tasks

Created by: Martin Sicho
On: 18-12-19, 13:14
"""
import time

from celery import shared_task

from .models import MolSet
from . import initializers

@shared_task(name='populateMolSet')
def populateMolSet(molset_id, initializer_class):
    instance = MolSet.objects.get(pk=molset_id)

    initializer_class = getattr(initializers, initializer_class)
    initializer = initializer_class(instance)
    initializer.populateInstance()
    time.sleep(86400)
    return {
        "populationSize" : 0
        , "errors" : []
    }