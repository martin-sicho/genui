from django.db import models
from djcelery_model.models import TaskMixin

from polymorphic.models import PolymorphicModel

from commons.models import TaskShortcutsMixIn, PolymorphicTaskManager
from compounds.models import MolSet, ActivitySet
from projects.models import DataSet

class ModelFileFormat(models.Model):
    fileExtension = models.CharField(max_length=32, blank=False, unique=True)
    description = models.TextField(max_length=10000, blank=True)

class Algorithm(models.Model):
    name = models.CharField(blank=False, max_length=128, unique=True)
    fileFormats = models.ManyToManyField(ModelFileFormat)

class ModelParameter(models.Model):
    STRING = 'string'
    BOOL = 'bool'
    INTEGER = 'integer'
    FLOAT = 'float'
    CONTENT_TYPES = [
       (STRING, 'String'),
       (BOOL, 'Logical'),
       (INTEGER, 'Integer'),
       (FLOAT, 'Float'),
    ]

    name = models.CharField(max_length=128, blank=False)
    algorithm = models.ForeignKey(Algorithm, on_delete=models.CASCADE, null=False, related_name='parameters')
    contentType = models.CharField(max_length=32, choices=CONTENT_TYPES, default=STRING)

    class Meta:
        unique_together = ('name', 'algorithm')

class ModelParameterValue(PolymorphicModel):
    parameter = models.ForeignKey(ModelParameter, on_delete=models.CASCADE, null=False)

class ModelParameterStr(ModelParameterValue):
    value = models.CharField(max_length=1024)

class ModelParameterBool(ModelParameterValue):
    value = models.BooleanField(null=False)

class ModelParameterInt(ModelParameterValue):
    value = models.IntegerField(null=False)

class ModelParameterFloat(ModelParameterValue):
    value = models.FloatField(null=False)

class ModelPerformanceMetric(models.Model):
    name = models.CharField(unique=True, blank=False, max_length=128)
    description = models.TextField(max_length=10000, blank=True)

class DescriptorGroup(models.Model):
    name = models.CharField(max_length=128, blank=False, unique=True)

class TrainingStrategy(PolymorphicModel):
    CLASSIFICATION = 'classification'
    REGRESSION = 'regression'
    MODES = [
       (CLASSIFICATION, 'Classification'),
       (REGRESSION, 'Regression'),
    ]

    algorithm = models.ForeignKey(Algorithm, on_delete=models.CASCADE, null=False)
    parameters = models.ManyToManyField(ModelParameterValue)
    fileFormat = models.ForeignKey(ModelFileFormat, on_delete=models.CASCADE, null=False)
    mode = models.CharField(choices=MODES, blank=False, max_length=32)

class QSARTrainingStrategy(TrainingStrategy):
    descriptors = models.ManyToManyField(DescriptorGroup)
    activityThreshold = models.FloatField(null=True)

class ValidationStrategy(PolymorphicModel):
    metrics = models.ManyToManyField(ModelPerformanceMetric)

class CV(ValidationStrategy):
    cvFolds = models.IntegerField(blank=False)

    class Meta:
        abstract = True

class ValidationSet(ValidationStrategy):
    validSetSize = models.FloatField(blank=False)

    class Meta:
        abstract = True

class BasicValidationStrategy(ValidationSet, CV):
    pass

class Model(TaskShortcutsMixIn, TaskMixin, DataSet):
    objects = PolymorphicTaskManager()

    trainingStrategy = models.ForeignKey(TrainingStrategy, null=False, on_delete=models.CASCADE)
    validationStrategy = models.ForeignKey(ValidationStrategy, null=False, on_delete=models.CASCADE)

class ModelPerformance(PolymorphicModel):
    metric = models.ForeignKey(ModelPerformanceMetric, null=False, on_delete=models.CASCADE)
    value = models.FloatField(blank=False)
    model = models.ForeignKey(Model, null=False, on_delete=models.CASCADE, related_name="performance")

class ModelPerformanceCV(ModelPerformance):
    fold = models.IntegerField(blank=False)

class ModelActivitySet(ActivitySet):
    pass

class QSARModel(Model):
    molset = models.ForeignKey(MolSet, null=False, on_delete=models.CASCADE, related_name="models")
    activities = models.ForeignKey(ModelActivitySet, null=True, on_delete=models.CASCADE, related_name="model")