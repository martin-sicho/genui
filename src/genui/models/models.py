import os

from django.db import models
import uuid

# Create your models here.
from djcelery_model.models import TaskMixin
from polymorphic.models import PolymorphicModel

from genui.utils.models import NON_POLYMORPHIC_CASCADE, OverwriteStorage
from genui.utils.extensions.tasks.models import TaskShortcutsMixIn, PolymorphicTaskManager
from genui.projects.models import DataSet


class AlgorithmMode(models.Model):
    name = models.CharField(unique=True, blank=False, max_length=32)

    def __str__(self):
        return '%s object (%s)' % (self.__class__.__name__, self.name)


class ModelFileFormat(models.Model):
    fileExtension = models.CharField(max_length=32, blank=False, unique=True)
    description = models.TextField(max_length=10000, blank=True)

class ImportableModelComponent(models.Model):

    corePackage = models.CharField(blank=False, null=False, default='genui.models.genuimodels', max_length=1024)

    class Meta:
        abstract = True


class Algorithm(ImportableModelComponent):
    name = models.CharField(blank=False, max_length=128, unique=True)
    fileFormats = models.ManyToManyField(ModelFileFormat)
    validModes = models.ManyToManyField(AlgorithmMode)

    def __str__(self):
        return '%s object (%s)' % (self.__class__.__name__, self.name)


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
    defaultValue = models.ForeignKey("ModelParameterValue", on_delete=models.SET_NULL, null=True)

    class Meta:
        unique_together = ('name', 'algorithm')

    def __str__(self):
        return '%s object (%s)' % (self.__class__.__name__, self.name)

class ModelBuilder(ImportableModelComponent):
    name = models.CharField(max_length=128, blank=False, unique=True)

    def __str__(self):
        return '%s object (%s)' % (self.__class__.__name__, self.name)

class ModelFile(models.Model):
    MAIN = "main"
    AUXILIARY = "aux"
    KINDS = [
       (MAIN, 'Main'),
       (AUXILIARY, 'Auxiliary'),
    ]

    class Rejected(Exception):

        def __init__(self, msg):
            super().__init__(msg)

    class InvalidFileFormatError(Exception):

        def __init__(self, msg):
            super().__init__(msg)

    modelInstance = models.ForeignKey("Model", null=False, related_name="files", on_delete=models.CASCADE)
    kind = models.CharField(max_length=32, choices=KINDS, null=False, default=AUXILIARY)
    note = models.CharField(max_length=128, blank=True)
    format = models.ForeignKey(ModelFileFormat, null=True, on_delete=models.CASCADE)
    file = models.FileField(null=True, upload_to='models/', storage=OverwriteStorage()) # TODO: add custom logic to save in a directory specific to the project where the model is

    @property
    def path(self):
        return self.file.path

    @staticmethod
    def generateMainFileName(model, fileFormat):
        return f"{model.trainingStrategy.algorithm.name}{model.id}_project{model.project.id}_{uuid.uuid4().hex}_main{fileFormat.fileExtension}"

    @staticmethod
    def generateAuxFileName(model, fileFormat):
        return f"{model.trainingStrategy.algorithm.name}{model.id}_project{model.project.id}_{uuid.uuid4().hex}_aux{fileFormat.fileExtension}"

    @staticmethod
    def create(model, name, file_, kind=AUXILIARY, note=None):
        if not note:
            note = ''
        if kind == ModelFile.MAIN and model.modelFile:
            algorithm = model.trainingStrategy.algorithm
            file_format = None
            for format_ in algorithm.fileFormats.all():
                if name.endswith(format_.fileExtension):
                    file_format = format_
                    break
            if not file_format:
                raise ModelFile.InvalidFileFormatError(f"The extension for file '{name}' of the submitted file did not match any of the known formats for algorithm: ({algorithm.name}).")

            if model.modelFile.format.fileExtension == file_format.fileExtension:
                model.modelFile.file.save(os.path.basename(model.modelFile.path), file_)
            else:
                model.modelFile.delete()
                ModelFile.objects.create(
                    model=model,
                    kind=ModelFile.MAIN,
                    format=file_format,
                    note=note,
                    file=file_
                )

            return model.modelFile
        else:
            file_format = None
            for format_ in ModelFileFormat.objects.all():
                if name.endswith(format_.fileExtension):
                    file_format = format_
                    break
            if kind == ModelFile.MAIN:
                if not file_format:
                    raise ModelFile.InvalidFileFormatError(f"The extension for file '{name}' of the submitted file did not match any of the known formats for model: {model.name}.")
                ret = ModelFile.objects.create(
                    modelInstance=model,
                    kind=ModelFile.MAIN,
                    format=file_format,
                    note=note
                )
                ret.file.save(ret.generateMainFileName(model, file_format), file_)
            else:
                ret = ModelFile.objects.create(
                    modelInstance=model,
                    kind=kind,
                    format=file_format if file_format else ModelFileFormat.objects.get_or_create(
                        fileExtension='.' + name.split('.')[-1]
                    )[0],
                    note=note
                )
                ret.file.save(ret.generateAuxFileName(model, ret.format), file_)

            return ret

class Model(TaskShortcutsMixIn, TaskMixin, DataSet):
    objects = PolymorphicTaskManager()
    builder = models.ForeignKey(ModelBuilder, on_delete=models.CASCADE, null=False)

    def __str__(self):
        return '%s object (%s)' % (self.__class__.__name__, self.name)

    @property
    def modelFile(self):
        # TODO: exception when more than one main file found
        main = self.files.filter(kind=ModelFile.MAIN)
        if main:
            return main.get()
        else:
            return None

    def onFileSave(self, saved : ModelFile):
        """
        This will be called when a file is being
        saved to this model instance. You can throw
        the ModelFile.Rejected exception if the file
        is invalid.

        :param saved:
        :return:
        """

        pass

    # @modelFile.setter
    # def modelFile(self, val):
    #     main = self.files.filter(kind=ModelFile.MAIN)
    #     if main:
    #         main.delete()
    #     val.kind = ModelFile.MAIN
    #     val.save()
    #     self.files.add(val)
    #     self.save()

    @property
    def trainingStrategy(self):
        count = self.trainingStrategies.count()
        if count == 1:
            return self.trainingStrategies.get()
        elif count == 0:
            return None
        else:
            raise Exception("Training strategy returned more than one value. This indicates an integrity error in the database!")

    # @property
    # def validationStrategy(self):
    #     count = self.validationStrategies.count()
    #     if count == 1:
    #         return self.validationStrategies.get()
    #     elif count == 0:
    #         return None
    #     else:
    #         raise Exception("Validation strategy returned more than one value. This indicates an integrity error in the database!")


class TrainingStrategy(PolymorphicModel):
    algorithm = models.ForeignKey(Algorithm, on_delete=models.CASCADE, null=False)
    mode = models.ForeignKey(AlgorithmMode, on_delete=models.CASCADE, null=False)
    modelInstance = models.ForeignKey(Model, null=False, on_delete=models.CASCADE, related_name="trainingStrategies")

    def processMetaData(self, metadata):
        pass


class ModelParameterValue(PolymorphicModel):
    parameter = models.ForeignKey(ModelParameter, on_delete=models.CASCADE, null=False)
    strategy = models.ForeignKey(TrainingStrategy, on_delete=NON_POLYMORPHIC_CASCADE, null=True, related_name='parameters')

    @staticmethod
    def parseValue(val):
        return str(val)


class ModelParameterStr(ModelParameterValue):
    value = models.CharField(max_length=1024)


class ModelParameterBool(ModelParameterValue):
    value = models.BooleanField(null=False)

    @staticmethod
    def parseValue(val):
        return bool(val)


class ModelParameterInt(ModelParameterValue):
    value = models.IntegerField(null=False)

    @staticmethod
    def parseValue(val):
        return int(val)


class ModelParameterFloat(ModelParameterValue):
    value = models.FloatField(null=False)

    @staticmethod
    def parseValue(val):
        return float(val)

PARAM_VALUE_CTYPE_TO_MODEL_MAP = {
    ModelParameter.STRING : ModelParameterStr,
    ModelParameter.INTEGER : ModelParameterInt,
    ModelParameter.FLOAT : ModelParameterFloat,
    ModelParameter.BOOL : ModelParameterBool
}


class ModelPerformanceMetric(ImportableModelComponent):
    name = models.CharField(unique=True, blank=False, max_length=128)
    validModes = models.ManyToManyField(AlgorithmMode, related_name='metrics')
    validAlgorithms = models.ManyToManyField(Algorithm, related_name='metrics')
    description = models.TextField(max_length=10000, blank=True)

    def __str__(self):
        return '%s object (%s)' % (self.__class__.__name__, self.name)


class ValidationStrategy(PolymorphicModel):
    metrics = models.ManyToManyField(ModelPerformanceMetric)
    trainingStrategy = models.ForeignKey(TrainingStrategy, null=False, on_delete=models.CASCADE, related_name='validationStrategies')
    # CHANGE: ValidationStrategy now linked to TrainingStrategy instead of Model.
    # This allows for more flexible configuration of validation strategies for different training approaches.

class CV(ValidationStrategy):
    cvFolds = models.IntegerField(blank=False)

    class Meta:
        abstract = True


class ValidationSet(ValidationStrategy):
    validSetSize = models.FloatField(blank=False)

    class Meta:
        abstract = True

# ako je ta basic si vytvorim validation a cv , pro kazdu abstraktnu si vytvorim 
# prislušnú konkretnu pre testovanie
class BasicValidationStrategy(ValidationSet, CV):
    pass


class ModelPerformance(PolymorphicModel):
    metric = models.ForeignKey(ModelPerformanceMetric, null=False, on_delete=models.CASCADE)
    value = models.FloatField(blank=False)
    model = models.ForeignKey(Model, null=False, on_delete=NON_POLYMORPHIC_CASCADE, related_name="performance")


class ModelPerformanceCV(ModelPerformance):
    fold = models.IntegerField(blank=False)

class ModelPerfomanceNN(ModelPerformance):
    epoch = models.IntegerField(null=False, blank=False)
    step = models.IntegerField(null=False, blank=False)

class ROCCurvePoint(ModelPerformance):

    fpr = models.FloatField(blank=False)
    auc = models.ForeignKey(ModelPerformance, null=False, on_delete=NON_POLYMORPHIC_CASCADE, related_name="points")

    @property
    def tpr(self):
        return self.value