from django.db import models
from djcelery_model.models import TaskMixin

from compounds.models import MolSet, ActivitySet
from modelling.models import ModelParameter, Model, TrainingStrategy


class DescriptorGroup(models.Model):
    name = models.CharField(max_length=128, blank=False, unique=True)

    def __str__(self):
        return '%s object (%s)' % (self.__class__.__name__, self.name)

class QSARTrainingStrategy(TrainingStrategy):
    descriptors = models.ManyToManyField(DescriptorGroup)
    activityThreshold = models.FloatField(null=True)


class QSARModel(Model):
    molset = models.ForeignKey(MolSet, null=True, on_delete=models.CASCADE, related_name="models")

class ModelActivitySet(ActivitySet):
    model = models.ForeignKey(QSARModel, null=False, on_delete=models.CASCADE, related_name="predictions")