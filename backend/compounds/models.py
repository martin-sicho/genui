from django_rdkit import models
from django_celery_results.models import TaskResult
from djcelery_model.models import TaskMixin, TaskManager
from polymorphic.managers import PolymorphicManager
from polymorphic.models import PolymorphicModel
from projects.models import DataSet

# Create your models here.

class PolymorphicTaskManager(PolymorphicManager, TaskManager):
    pass

class TaskShortcutsMixIn:

    def getTasksAsDict(self, started_only=False):
        if started_only:
            tasks = self.tasks.started()
        else:
            tasks = self.tasks.all()

        grouped_tasks = dict()
        for task in tasks:
            task_id = task.task_id
            try:
                result = TaskResult.objects.get(task_id=task_id)
                task_name = result.task_name
            except TaskResult.DoesNotExist:
                print(f"Task {task_id} not found in the database. Skipping...")
                continue
            if not task_name:
                task_name = 'UnknownTask'
            if task_name not in grouped_tasks:
                grouped_tasks[task_name] = []
            grouped_tasks[task_name].append(result)

        data = dict()
        for key in grouped_tasks:
            data[key] = [{"task_id" : x.task_id, "status" : x.status} for x in grouped_tasks[key]]
        return data

class MolSet(TaskShortcutsMixIn, TaskMixin, DataSet):
    objects = PolymorphicTaskManager()

class ActivitySet(TaskShortcutsMixIn, TaskMixin, DataSet):
    objects = PolymorphicTaskManager()

    molecules = models.ForeignKey(MolSet, blank=False, null=True, on_delete=models.CASCADE, related_name="activities") # FIXME: it probably makes more sense to make this field non-nullable

class Molecule(PolymorphicModel):
    molObject = models.MolField()
    canonicalSMILES = models.CharField(max_length=65536)
    inchiKey = models.CharField(max_length=65536, unique=True)
    providers = models.ManyToManyField(MolSet, blank=False, related_name='molecules')

class ChEMBLAssay(models.Model):
    assayID = models.CharField(max_length=32, unique=True, blank=False)

class ChEMBLTarget(models.Model):
    targetID = models.CharField(max_length=32, unique=True, blank=False)

class ChEMBLMolecule(Molecule):
    chemblID = models.CharField(max_length=32, unique=True, blank=False, null=False)

class ChEMBLActivities(ActivitySet):
    pass

class ChEMBLCompounds(MolSet):
    targets = models.ManyToManyField(ChEMBLTarget, blank=False)

class ActivityUnits(models.Model):
    value = models.CharField(blank=False, max_length=8, unique=True)

class Activity(models.Model):
    value = models.FloatField(blank=False)
    units = models.ForeignKey(ActivityUnits, on_delete=models.CASCADE, null=True)
    source = models.ForeignKey(ActivitySet, on_delete=models.CASCADE, blank=False)
    molecule = models.ForeignKey(Molecule, on_delete=models.CASCADE, blank=False)

class ChEMBLActivity(Activity):
    type = models.CharField(blank=False, max_length=128)
    relation = models.CharField(blank=False, max_length=128)
    assay = models.ForeignKey(ChEMBLAssay, on_delete=models.CASCADE, null=False, blank=False)
    target = models.ForeignKey(ChEMBLTarget, on_delete=models.CASCADE, null=False, blank=False)
    comment = models.CharField(blank=True, max_length=128, null=True)
