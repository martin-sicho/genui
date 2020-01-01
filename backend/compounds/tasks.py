"""
tasks

Created by: Martin Sicho
On: 18-12-19, 13:14
"""
from celery import shared_task
from celery_progress.backend import ProgressRecorder

from .models import MolSet
from .initializers import populate_molset

@shared_task(name='CreateCompoundSet', bind=True)
def populateMolSet(self, molset_id, initializer_class, initializer_kwargs=None):
    if not initializer_kwargs:
        initializer_kwargs = dict()
    instance = MolSet.objects.get(pk=molset_id)
    initializer_class = getattr(populate_molset, initializer_class)
    initializer = initializer_class(instance, ProgressRecorder(self), **initializer_kwargs)
    count = initializer.populateInstance()
    return {
        "populationSize" : count
        , "errors" : [repr(x) for x in initializer.errors]
    }