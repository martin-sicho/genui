"""
tasks

Created by: Martin Sicho
On: 18-12-19, 13:14
"""
from celery import shared_task

from genui.utils.inspection import getObjectAndModuleFromFullName
from genui.utils.extensions.tasks.progress import ProgressRecorder
from genui.compounds.models import MolSet, MolSetExport


def populate(self, molset_id, initializer_class, initializer_kwargs=None):
    if not initializer_kwargs:
        initializer_kwargs = dict()
    instance = MolSet.objects.get(pk=molset_id)
    initializer_class, initializer_module = getObjectAndModuleFromFullName(initializer_class)
    initializer = initializer_class(instance, ProgressRecorder(self), **initializer_kwargs)
    count = initializer.populateInstance()
    return {
        "populationSize" : count
        , "errors" : [repr(x) for x in initializer.errors]
    }

@shared_task(name='CreateCompoundSet', bind=True)
def populateMolSet(self, molset_id, initializer_class, initializer_kwargs=None):
    return populate(self, molset_id, initializer_class, initializer_kwargs)

@shared_task(name='CreateCompoundSetGPU', bind=True, queue='gpu')
def populateMolSetGPU(self, molset_id, initializer_class, initializer_kwargs=None):
    return populate(self, molset_id, initializer_class, initializer_kwargs)

@shared_task(name='UpdateCompoundSet', bind=True)
def updateMolSet(self, molset_id, updater_class, updater_kwargs=None):
    if not updater_kwargs:
        updater_kwargs = dict()
    instance = MolSet.objects.get(pk=molset_id)
    updater_class, updater_module = getObjectAndModuleFromFullName(updater_class)
    updater = updater_class(instance, ProgressRecorder(self), **updater_kwargs)
    count = updater.updateInstance()
    return {
        "populationSize" : count
        , "errors" : [repr(x) for x in updater.errors]
    }

@shared_task(name='CreateExport', bind=True)
def createExport(self, export_id):
    instance = MolSetExport.objects.get(pk=export_id)
    exporter_class, exporter_module = getObjectAndModuleFromFullName(instance.exporter.classPath)
    exporter = exporter_class(instance, ProgressRecorder(self))
    file = exporter.saveFile()
    return {
        "file" : file.file.name
        , "errors" : [repr(x) for x in exporter.errors]
    }