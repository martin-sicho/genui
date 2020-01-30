"""
bases

Created by: Martin Sicho
On: 24-01-20, 15:03
"""
import traceback
from abc import ABC, abstractmethod

import joblib
from pandas import DataFrame, Series

import uuid
from django.core.files.base import ContentFile

from modelling import models
from commons.helpers import findClassInModule


class Algorithm(ABC):
    name = None
    CLASSIFICATION = 'classification'
    REGRESSION = 'regression'

    @classmethod
    def attachModesToModel(cls, model, modes):
        model.validModes.clear()
        for mode in modes:
            mode = models.AlgorithmMode.objects.get_or_create(name=mode)[0]
            model.validModes.add(mode)
        model.save()
        return modes

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
    def getDjangoModel(cls) -> models.Algorithm:
        if not cls.name:
            raise Exception('You have to specify a name for the algorithm in its class "name" property')
        ret = models.Algorithm.objects.get_or_create(
            name=cls.name
        )[0]

        cls.attachModesToModel(ret, cls.getModes()) # TODO: this should use the same pattern as the file formats method
        cls.getFileFormats(attach_to=ret)
        return ret

    @staticmethod
    @abstractmethod
    def getParams():
        pass

    def __init__(self, builder, callback=None):
        self.builder = builder
        self.trainingInfo = builder.training
        self.validationInfo = builder.validation
        self.params = {x.parameter.name : x.value for x in self.trainingInfo.parameters.all()}
        self.mode = self.trainingInfo.mode
        self.callback = callback
        self._model = None

    def getSerializer(self):
        return lambda filename : joblib.dump(
            self._model
            , filename
        )

    def serialize(self, filename):
        self.getSerializer()(filename)

    def deserialize(self, filename):
        self._model = joblib.load(filename)

    @property
    @abstractmethod
    def model(self):
        pass

    @abstractmethod
    def fit(self, X : DataFrame, y : Series):
        pass

    @abstractmethod
    def predict(self, X : DataFrame):
        pass


class ValidationMetric(ABC):
    name = None
    description = None

    def __init__(self, builder):
        self.builder = builder

    @classmethod
    def getDjangoModel(cls):
        if not cls.name:
            raise Exception('You have to specify a name for the validation metric in its class "name" property')
        ret = models.ModelPerformanceMetric.objects.get_or_create(
            name=cls.name
        )[0]
        if cls.description:
            ret.description = cls.description
            ret.save()
        return ret

    @abstractmethod
    def __call__(self, true_vals : Series, predicted_vals : Series):
        pass

    def save(
            self,
            true_vals : Series,
            predicted_vals : Series,
            perfClass=models.ModelPerformance,
            **kwargs
    ):
        perfClass.objects.create(
                    metric=models.ModelPerformanceMetric.objects.get(name=self.name),
                    value=self(true_vals, predicted_vals),
                    model=self.builder.instance,
                    **kwargs
                )

class ModelBuilder(ABC):

    @staticmethod
    def findAlgorithmClass(name):
        from . import algorithms
        return findClassInModule(Algorithm, algorithms, "name", name)

    @staticmethod
    def findMetricClass(name):
        from . import metrics
        return findClassInModule(ValidationMetric, metrics, "name", name)

    def __init__(
            self,
            instance : models.Model,
            progress = None,
            onFit = None
    ):
        self.instance = instance

        self.training = self.instance.trainingStrategy
        self.algorithmClass = self.findAlgorithmClass(self.training.algorithm.name)
        self.onFit = onFit

        self.validation = self.instance.validationStrategy
        self.metricClasses = [self.findMetricClass(x.name) for x in self.validation.metrics.all()]

        self.progress = progress
        self.progressStages = []
        self.currentProgress = 0
        self.errors = []

        self._model = self.algorithmClass(self, self.onFit).deserialize(self.instance.modelFile.path) if self.instance.modelFile else None


    @property
    def model(self) -> Algorithm:
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

    def recordProgress(self):
        if self.currentProgress < len(self.progressStages):
            if self.progress:
                self.progress.set_progress(
                    self.currentProgress
                    , len(self.progressStages)
                    , description=self.progressStages[self.currentProgress]
                )
            print(self.progressStages[self.currentProgress])
        else:
            self.errors.append(Exception("Incorrect progress count detected."))
        self.currentProgress += 1
        print(f"{self.currentProgress}/{len(self.progressStages)}")

    def build(self) -> models.Model:
        self.model = self.algorithmClass(self, self.onFit)
        self.model.fit(self.getX(), self.getY())
        self.saveFile()
        return self.instance

    def saveFile(self):
        name = f"{self.algorithmClass.name}{self.instance.id}_project{self.instance.project.id}_{uuid.uuid1()}"
        if not self.instance.modelFile:
            model_format = self.training.algorithm.fileFormats.all()[0] # FIXME: this should be changed once we expose the file formats in the training strategy
            self.instance.modelFile.save(name + model_format.fileExtension, ContentFile('Dummy file for {0}'.format(name)))
        path = self.instance.modelFile.path
        self.model.serialize(path)

    def predict(self, X : DataFrame):
        if self.model:
            self.model.predict(X)
        else:
            raise Exception("The model is not trained or loaded. Invalid call to 'predict'.") # TODO: throw a more specific exception

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
        for metric_class in self.metricClasses:
            try:
                metric_class(self).save(y_validated, y_predicted, perfClass, *args, **kwargs)
            except Exception as exp:
                # TODO: add special exception
                print("Failed to obtain values for metric: ", metric_class.name)
                self.errors.append(exp)
                traceback.print_exc()
                continue

class ModelBuilderWithValidation(ValidationMixIn, ModelBuilder, ABC):
    pass