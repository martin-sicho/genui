from django.db import models
from django_celery_results.models import TaskResult
from djcelery_model.models import TaskMixin, TaskManager
from polymorphic.managers import PolymorphicManager
from polymorphic.models import PolymorphicModel
from projects.models import DataSet

# Create your models here.

class PolymorphicTaskManager(PolymorphicManager, TaskManager):
    pass

class MolSet(TaskMixin, DataSet):
    objects = PolymorphicTaskManager()

    def getTasksAsDict(self, started_only=False):
        if started_only:
            tasks = self.tasks.started()
        else:
            tasks = self.tasks.all()

        grouped_tasks = dict()
        for task in tasks:
            task_id = task.task_id
            result = TaskResult.objects.get(task_id=task_id)
            task_name = result.task_name
            if task_name not in grouped_tasks:
                grouped_tasks[task_name] = []
            grouped_tasks[task_name].append(result)

        data = dict()
        for key in grouped_tasks:
            data[key] = [{"task_id" : x.task_id, "status" : x.status} for x in grouped_tasks[key]]
        return data

class ActivitySet(DataSet):
    pass

class Molecule(PolymorphicModel):

    canonicalSMILES = models.CharField(max_length=65536)
    inchiKey = models.CharField(max_length=65536, unique=True)
    providers = models.ManyToManyField(MolSet, blank=False, related_name='molecules')

class ChEMBLAssay(models.Model):
    assayID = models.CharField(max_length=32, unique=True)

class ChEMBLTarget(models.Model):
    targetID = models.CharField(max_length=32, unique=True)

class ChEMBLMolecule(Molecule):
    chemblID = models.CharField(max_length=32, unique=True, blank=False, null=False)
    assays = models.ManyToManyField(ChEMBLAssay, blank=True)
    targets = models.ManyToManyField(ChEMBLTarget, blank=True)

class ChEMBLCompounds(MolSet):
    assays = models.ManyToManyField(ChEMBLAssay, blank=True)
    targets = models.ManyToManyField(ChEMBLTarget, blank=True)

class ChEMBLActivities(ActivitySet):
    assays = models.ManyToManyField(ChEMBLAssay, blank=False)
    targets = models.ManyToManyField(ChEMBLTarget, blank=False)

class ActivityUnit(models.Model):
    value = models.CharField(blank=False, max_length=8, unique=True)

class Activity(models.Model):

    value = models.FloatField(blank=False)
    units = models.ForeignKey(ActivityUnit, on_delete=models.CASCADE, blank=False)
    source = models.ForeignKey(ActivitySet, on_delete=models.CASCADE, blank=False)
    molecule = models.ForeignKey(Molecule, on_delete=models.CASCADE, blank=False)

