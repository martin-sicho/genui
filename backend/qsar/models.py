from django.db import models
from djcelery_model.models import TaskMixin

from polymorphic.models import PolymorphicModel

from commons.models import TaskShortcutsMixIn, PolymorphicTaskManager
from projects.models import DataSet

class ModelFileFormat(models.Model):
    fileExtension = models.CharField(max_length=8, blank=False, unique=True)
    description = models.TextField(max_length=10000, blank=True)

class Algorithm(models.Model):
    name = models.CharField(blank=False, max_length=128, unique=True)
    fileFormats = models.ManyToManyField(ModelFileFormat)

class ModelParameters(PolymorphicModel):
    algorithm = models.ForeignKey(Algorithm, null=False, on_delete=models.CASCADE)

class RandomForestParams(ModelParameters):
    nTrees = models.IntegerField(blank=False)
    # TODO: complete this

class ModelPerformanceMetric(models.Model):
    name = models.CharField(unique=True, blank=False, max_length=128)
    description = models.TextField(max_length=10000, blank=True)

class TrainingStrategy(models.Model):
    name = models.CharField(max_length=128, unique=True)
    supportedAlgorithms = models.ManyToManyField(Algorithm, related_name="trainingStrategies")
    supportedMetrics = models.ManyToManyField(ModelPerformanceMetric, related_name="metrics")

class TrainingStrategyImpl(PolymorphicModel):
    strategy = models.ForeignKey(TrainingStrategy, on_delete=models.CASCADE, null=False)

class CVStrategy(TrainingStrategy):
    cvFolds = models.IntegerField(blank=False)

    class Meta:
        abstract = True

class ValidationSetStrategy(TrainingStrategy):
    validSetSize = models.FloatField(blank=False)

    class Meta:
        abstract = True

class BasicTrainingStrategy(ValidationSetStrategy, CVStrategy):
    pass

class Model(TaskShortcutsMixIn, TaskMixin, DataSet):
    objects = PolymorphicTaskManager()

    parameters = models.ForeignKey(ModelParameters, null=False, on_delete=models.CASCADE)
    training = models.ForeignKey(TrainingStrategyImpl, null=False, on_delete=models.CASCADE)

class ModelPerformance(PolymorphicModel):
    metric = models.ForeignKey(ModelPerformanceMetric, null=False, on_delete=models.CASCADE)
    value = models.FloatField(blank=False)
    model = models.ForeignKey(Model, null=False, on_delete=models.CASCADE, related_name="performance")

class ModelPerformanceCV(ModelPerformance):
    fold = models.IntegerField(blank=False)
    performance = models.ForeignKey(ModelPerformance, null=False, on_delete=models.CASCADE, related_name='performanceCV')