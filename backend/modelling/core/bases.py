"""
bases

Created by: Martin Sicho
On: 24-01-20, 15:03
"""
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

    @staticmethod
    def getModes():
        return [Algorithm.CLASSIFICATION, Algorithm.REGRESSION]

    @classmethod
    def getDjangoModel(cls) -> models.Algorithm:
        if not cls.name:
            raise Exception('You have to specify a name for the algorithm in its class "name" property')
        ret = models.Algorithm.objects.get_or_create(
            name=cls.name
        )[0]
        file_format = models.ModelFileFormat.objects.get_or_create(
            fileExtension=".joblib.gz",
            description="A compressed joblib file."
        )[0]
        ret.fileFormats.add(file_format)
        for mode in cls.getModes():
            mode = models.AlgorithmMode.objects.get_or_create(name=mode)[0]
            ret.validModes.add(mode)
        ret.save()
        return ret

    @staticmethod
    @abstractmethod
    def getParams():
        pass

    def __init__(self, builder, callback=None):
        self.trainingInfo = builder.training
        self.params = {x.parameter.name : x.value for x in self.trainingInfo.parameters.all()}
        self.mode = self.trainingInfo.mode
        self.fileFormat = self.trainingInfo.algorithm.fileFormats.all()[0]
        self.callback = callback
        self._model = None

    def getSerializer(self):
        return lambda filename : joblib.dump(
            self._model
            , filename + self.fileFormat.fileExtension
        )

    def serialize(self, filename):
        self.getSerializer()(filename)

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

        self._model = None


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

    def fit(self, callback=None) -> models.Model:
        # TODO: sanity check if length of X and y are the same
        self.model = self.algorithmClass(self, callback if callback else self.onFit)
        self.model.fit(self.getX(), self.getY())
        self.saveFile(self.model)
        return self.instance

    def saveFile(self, model : Algorithm):
        name = f"{self.algorithmClass.name}{self.instance.id}_project{self.instance.project.id}_{uuid.uuid1()}"
        extension = model.fileFormat.fileExtension
        self.instance.modelFile.save(name + extension, ContentFile('Dummy file for {0}'.format(name)))
        path = self.instance.modelFile.path.replace(extension, '')
        model.serialize(path)