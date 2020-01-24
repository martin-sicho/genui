from django.db import models
from djcelery_model.models import TaskMixin

from compounds.models import MolSet, ActivitySet
from modelling.models import ModelParameter, Model, TrainingStrategy, ModelParameterStr, ModelParameterBool, \
    ModelParameterInt, ModelParameterFloat

PARAM_VALUE_CTYPE_TO_MODEL_MAP = {
    ModelParameter.STRING : ModelParameterStr,
    ModelParameter.INTEGER : ModelParameterInt,
    ModelParameter.FLOAT : ModelParameterFloat,
    ModelParameter.BOOL : ModelParameterBool
}

class DescriptorGroup(models.Model):
    name = models.CharField(max_length=128, blank=False, unique=True)

class QSARTrainingStrategy(TrainingStrategy):
    descriptors = models.ManyToManyField(DescriptorGroup)
    activityThreshold = models.FloatField(null=True)


class QSARModel(Model):
    molset = models.ForeignKey(MolSet, null=False, on_delete=models.CASCADE, related_name="models")

class ModelActivitySet(ActivitySet):
    model = models.ForeignKey(QSARModel, null=False, on_delete=models.CASCADE, related_name="predictions")