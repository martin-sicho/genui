"""
bases

Created by: Martin Sicho
On: 24-01-20, 15:03
"""
from abc import ABC, abstractmethod

import joblib
from pandas import DataFrame, Series

import modelling.models


class Algorithm(ABC):
    name = None
    CLASSIFICATION = 'classification'
    REGRESSION = 'regression'

    @staticmethod
    def getModes():
        return [Algorithm.CLASSIFICATION, Algorithm.REGRESSION]

    @classmethod
    def getDjangoModel(cls) -> modelling.models.Algorithm:
        if not cls.name:
            raise Exception('You have to specify a name for the algorithm in its class "name" property')
        ret = modelling.models.Algorithm.objects.get_or_create(
            name=cls.name
        )[0]
        file_format = modelling.models.ModelFileFormat.objects.get_or_create(
            fileExtension=".joblib.gz",
            description="A compressed joblib file."
        )[0]
        ret.fileFormats.add(file_format)
        for mode in cls.getModes():
            mode = modelling.models.AlgorithmMode.objects.get_or_create(name=mode)[0]
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
        self.model = builder.instance

    @classmethod
    def getDjangoModel(cls):
        if not cls.name:
            raise Exception('You have to specify a name for the validation metric in its class "name" property')
        ret = modelling.models.ModelPerformanceMetric.objects.get_or_create(
            name=cls.name
        )[0]
        if cls.description:
            ret.description = cls.description
            ret.save()
        return ret

    @abstractmethod
    def __call__(self, true_vals : Series, predicted_vals : Series):
        pass