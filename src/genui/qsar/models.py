from django.db import models

from genui.compounds.models import MolSet, ActivitySet, Activity, ActivityTypes, ActivityUnits
from genui.models.models import Model, TrainingStrategy

from genui.models.models import ImportableModelComponent

class DescriptorGroup(ImportableModelComponent):
    name = models.CharField(max_length=128, blank=False, unique=True)

    def __str__(self):
        return '%s object (%s)' % (self.__class__.__name__, self.name)

class QSARTrainingStrategy(TrainingStrategy):
    descriptors = models.ManyToManyField(DescriptorGroup)
    activityThreshold = models.FloatField(null=True)
    activitySet = models.ForeignKey(ActivitySet, null=True, on_delete=models.CASCADE)
    activityType = models.ForeignKey(ActivityTypes, on_delete=models.CASCADE, null=True)

class QSARModel(Model):
    molset = models.ForeignKey(MolSet, null=True, on_delete=models.CASCADE, related_name="models")
    predictionsType = models.ForeignKey(ActivityTypes, on_delete=models.CASCADE, null=True)
    predictionsUnits = models.ForeignKey(ActivityUnits, on_delete=models.CASCADE, null=True)

class ModelActivitySet(ActivitySet):
    model = models.ForeignKey(QSARModel, null=False, on_delete=models.CASCADE, related_name="predictions")

class ModelActivity(Activity):
    pass