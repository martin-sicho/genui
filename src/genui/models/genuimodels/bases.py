"""
bases

Created by: Martin Sicho
On: 24-01-20, 15:03
"""
import traceback
import weakref
from abc import ABC, abstractmethod

import joblib
from django.core.exceptions import ImproperlyConfigured
from django.db import transaction
from pandas import DataFrame, Series

from django.core.files.base import ContentFile

from genui.models import models
from genui.utils.inspection import findSubclassByID, importFromPackage
from genui.models.models import ModelFile


class Algorithm(ABC):
    # TODO: use a metaclass to initialize these (the classmethods below have to be called explicitly at the moment)
    name = None
    parameters = {}
    CLASSIFICATION = 'classification'
    REGRESSION = 'regression'
    GENERATOR = 'generator'
    MAP = 'map'

    django_model = None
    django_parameters = []
    django_file_formats = []
    django_modes = []

    @classmethod
    def attachModesToModel(cls, model, modes):
        model.validModes.clear()
        for mode in modes:
            mode = models.AlgorithmMode.objects.get_or_create(name=mode)[0]
            model.validModes.add(mode)
        model.save()
        return model.validModes

    @classmethod
    def getModes(cls):
        return [cls.CLASSIFICATION, cls.REGRESSION]

    @staticmethod
    def attachToInstance(instance, items, field):
        field.clear()
        for item in items:
            field.add(item)
        instance.save()

    @classmethod
    def getFileFormats(cls, attach_to=None):
        formats = [models.ModelFileFormat.objects.get_or_create(
            fileExtension=".joblib.gz",
            description="A compressed joblib file."
        )[0]]

        if attach_to:
            cls.attachToInstance(attach_to, formats, attach_to.fileFormats)
        return formats

    @classmethod
    def getDjangoModel(cls, corePackage=None, update=False) -> models.Algorithm or None:
        # TODO: this should go to the init of the metaclass
        if not cls.name:
            print('This class has invalid name attribute. No django model can be provided for: ', cls.__name__)
            return

        ret, ret_created = models.Algorithm.objects.get_or_create(name=cls.name)

        # just return if we are not setting up a new instance
        if not ret_created and not update:
            return ret

        if corePackage:
            ret.corePackage = corePackage
            ret.save()

        cls.django_modes = cls.attachModesToModel(ret, cls.getModes()) # TODO: this should use the same pattern as the file formats method
        cls.django_file_formats = cls.getFileFormats(attach_to=ret)
        cls.django_model = ret
        cls.django_parameters = cls.getParams()
        return ret

    @classmethod
    def getParams(cls):
        ret = []
        for param_name in cls.parameters:
            param_type = cls.parameters[param_name]["type"]
            with transaction.atomic():
                param, created = models.ModelParameter.objects.get_or_create(
                    name=param_name,
                    contentType=param_type,
                    algorithm=cls.django_model
                )

                default_value = cls.parameters[param_name]["defaultValue"]
                if created:
                    print(f"Creating default value: {cls.__name__}.{param_name} = {default_value}")
                    param.defaultValue = models.PARAM_VALUE_CTYPE_TO_MODEL_MAP[param_type].objects.create(
                        parameter=param,
                        value=default_value
                    )
                    param.save()
                else:
                    if param.defaultValue.value != default_value:
                        print(f'Changing default value: {param_name}={default_value}')
                        param.defaultValue.value = default_value
                        param.defaultValue.save()

                ret.append(param)
        return ret

    def __init__(self, builder, callback=None):
        self._builder = weakref.ref(builder)
        self.instance = self.builder.instance
        self.trainingInfo = self.builder.training
        self.validationInfo = self.builder.validation
        self.params = {x.parameter.name : x.value for x in self.trainingInfo.parameters.all()}
        self.mode = self.trainingInfo.mode
        self.callback = callback
        self._model = None

    @property
    def builder(self):
        ret = self._builder()
        if ret:
            return ret
        else:
            raise LookupError("Builder was destroyed before being referenced!")

    def getSerializer(self):
        return lambda filename : joblib.dump(
            self.model
            , filename
        )

    def serialize(self, filename):
        self.getSerializer()(filename)

    def getDeserializer(self):
        return lambda filename : joblib.load(filename)

    def deserialize(self, filename):
        self._model = self.getDeserializer()(filename)
        return self

    @property
    @abstractmethod
    def model(self):
        pass

    @abstractmethod
    def fit(self, X : DataFrame, y : Series):
        pass

    @abstractmethod
    def predict(self, X : DataFrame) -> Series:
        pass


class ValidationMetric(ABC):
    name = None
    description = None

    def __init__(self, builder):
        self.builder = builder

    @classmethod
    def getDjangoModel(cls, corePackage=None, update=False):
        if not cls.name:
            raise Exception('You have to specify a name for the validation metric in its class "name" property')

        ret, ret_created = models.ModelPerformanceMetric.objects.get_or_create(
            name=cls.name
        )

        # just return if we are not setting up a new instance
        if not ret_created and not update:
            return ret

        # just return if we are not creating a new instance
        if corePackage:
            ret.corePackage = corePackage
            ret.save()
        if hasattr(cls, 'description'):
            ret.description = cls.description
            ret.save()

        ret.validModes.clear()
        if hasattr(cls, 'modes'):
            for mode in cls.modes:
                mode = models.AlgorithmMode.objects.get_or_create(
                    name=mode
                )[0]
                ret.validModes.add(mode)

        ret.validAlgorithms.clear()
        if hasattr(cls, 'algorithms'):
            for alg in cls.algorithms:
                alg = models.Algorithm.objects.get_or_create(
                    name=alg.name
                )[0]
                ret.validAlgorithms.add(alg)
        return ret

    @abstractmethod
    def __call__(self, true_vals : Series, predicted_vals : Series):
        pass

    @staticmethod
    def probasToClasses(probas):
        return [1 if x >= 0.5 else 0 for x in probas]

    def save(
            self,
            true_vals : Series,
            predicted_vals : Series,
            perfClass=models.ModelPerformance,
            **kwargs
    ):
        return perfClass.objects.create(
                    metric=models.ModelPerformanceMetric.objects.get(name=self.name),
                    value=self(true_vals, predicted_vals),
                    model=self.builder.instance,
                    **kwargs
                )

class ModelBuilder(ABC):

    @classmethod
    def getDjangoModel(cls, corePackage=None, update=False):
        ret, ret_created = models.ModelBuilder.objects.get_or_create(
            name=cls.__name__
        )

        # just return if we are not setting up a new instance
        if not ret_created and not update:
            return ret

        if corePackage:
            ret.corePackage = corePackage
            ret.save()
        return ret

    def findAlgorithmClass(self, name, corePackage=None):
        if not corePackage:
            corePackage = self.corePackage
        return findSubclassByID(
            Algorithm,
            importFromPackage(corePackage, "algorithms"),
            "name",
            name
        )

    def findMetricClass(self, name, corePackage=None):
        if not corePackage:
            corePackage = self.corePackage
        return findSubclassByID(
            ValidationMetric
            , importFromPackage(corePackage, "metrics")
            , "name"
            , name
        )

    def __init__(
            self,
            instance : models.Model,
            progress = None,
            onFit = None
    ):
        self.instance = instance

        self.training = self.instance.trainingStrategy
        self.algorithmPackage = self.training.algorithm.corePackage
        self.algorithmClass = self.findAlgorithmClass(
            self.training.algorithm.name,
            self.algorithmPackage if self.algorithmPackage else None
        )
        self.onFit = onFit

        self.validation = self.instance.validationStrategy
        self.metricClasses = [self.findMetricClass(x.name, x.corePackage) for x in self.validation.metrics.all()] if self.validation else []

        self.progress = progress
        self.errors = []

        self._model = None

    @property
    def corePackage(self):
        from .. import genuimodels
        return genuimodels

    @property
    def model(self) -> Algorithm:
        if self._model is None:
            self._model = self.algorithmClass(self, self.onFit)
            if self.instance.modelFile:
                self._model.deserialize(self.instance.modelFile.path)
        return self._model

    @model.setter
    def model(self, val):
        self._model = val

    @abstractmethod
    def getY(self) -> Series:
        pass

    @abstractmethod
    def getX(self) -> DataFrame:
        pass

    def build(self) -> models.Model:
        self.model.fit(self.getX(), self.getY())
        self.saveFile()
        return self.instance

    def saveFile(self):
        if not self.instance.modelFile:
            model_format = self.training.algorithm.fileFormats.all()[0] # FIXME: this should be changed once we expose the file formats in the training strategy
            ModelFile.create(
                self.instance,
                f'main.{model_format.fileExtension}',
                ContentFile('placeholder'),
                kind=ModelFile.MAIN,
                note=f'{self.training.algorithm.name}_main'
            )
        path = self.instance.modelFile.path
        self.model.serialize(path)

class ProgressMixIn:

    def __init__(self, instance, progress, *args, **kwargs):
        super().__init__(instance, progress, *args, **kwargs)

        self.progress = progress
        self.progressStages = []
        self.currentProgress = 0
        self.errors = self.errors if hasattr(self, "errors") else []

    def recordProgress(self):
        if self.currentProgress < len(self.progressStages):
            if self.progress:
                self.progress.set_progress(
                    self.currentProgress+1
                    , len(self.progressStages)
                    , description=self.progressStages[self.currentProgress]
                )
            print(self.progressStages[self.currentProgress])
        else:
            self.errors.append(Exception("Incorrect progress count detected."))
        self.currentProgress += 1
        print(f"{self.currentProgress}/{len(self.progressStages)}")

class ValidationMixIn:

    def fitAndValidate(
            self,
            X_train : DataFrame,
            y_train : Series,
            X_validated : DataFrame,
            y_validated : Series,
            y_predicted=None,
            perfClass=models.ModelPerformance,
            *args,
            **kwargs
    ):
        if not y_predicted:
            model = self.algorithmClass(self)
            model.fit(X_train, y_train)
            y_predicted = model.predict(X_validated)
        self.validate(y_validated, y_predicted, perfClass, *args, **kwargs)


    def validate(
            self,
            y_validated,
            y_predicted,
            perfClass=models.ModelPerformance,
            *args,
            **kwargs):
        if not self.validation:
            raise ImproperlyConfigured(f"No validation strategy is set for model: {repr(self.instance)}")
        for metric_class in self.metricClasses:
            try:
                metric_class(self).save(y_validated, y_predicted, perfClass, *args, **kwargs)
            except Exception as exp:
                # TODO: add special exception
                print("Failed to obtain values for metric: ", metric_class.name)
                self.errors.append(exp)
                traceback.print_exc()
                continue

class PredictionMixIn:

    def predict(self, X : DataFrame = None) -> Series:
        if X is None:
            X = self.getX()
        # TODO: check if X is valid somehow
        if self.model:
            return self.model.predict(X)
        else:
            raise Exception("The model is not trained or loaded. Invalid call to 'predict'.") # TODO: throw a more specific exception


class CompleteBuilder(PredictionMixIn, ValidationMixIn, ProgressMixIn, ModelBuilder, ABC):
    pass