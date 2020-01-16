"""
algorithms

Created by: Martin Sicho
On: 14-01-20, 10:16
"""
from abc import ABC, abstractmethod

from commons.helpers import findClassInModule
from qsar import models
from pandas import DataFrame, Series


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
        for mode in Algorithm.getModes():
            mode = models.AlgorithmMode.objects.get_or_create(name=mode)[0]
            ret.validModes.add(mode)
        ret.save()
        return ret

    @staticmethod
    @abstractmethod
    def getParams():
        pass

    def __init__(self, training_info : models.TrainingStrategy, callback=None):
        self.trainingInfo = training_info
        self.params = {x.parameter.name : x.value for x in self.trainingInfo.parameters.all()}
        self.mode = self.trainingInfo.mode
        self.file = self.trainingInfo.fileFormat
        self.callback = callback

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

    @abstractmethod
    def __call__(self, true_vals : Series, predicted_vals : Series):
        pass

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

class DescriptorCalculator:
    group_name = None

    def __call__(self, smiles):
        pass

    @classmethod
    def getDjangoModel(cls) -> models.DescriptorGroup:
        if not cls.group_name:
            raise Exception('You have to specify a name for the descriptor group in its class "group_name" property')
        return models.DescriptorGroup.objects.get_or_create(name=cls.group_name)[0]


class QSARModelBuilder:
    @staticmethod
    def findAlgorithmClass(name):
        from . import algorithms
        return findClassInModule(Algorithm, algorithms, "name", name)

    @staticmethod
    def findMetricClass(name):
        from . import metrics
        return findClassInModule(ValidationMetric, metrics, "name", name)

    @staticmethod
    def findDescriptorClass(name):
        from . import descriptors
        return findClassInModule(DescriptorCalculator, descriptors, "group_name", name)

    def __init__(
            self,
            instance : models.QSARModel,
            progress = None,
            onFitCall=None
    ):
        self.instance = instance
        self.training = self.instance.trainingStrategy
        self.validation = self.instance.validationStrategy
        self.descriptorClasses = [self.findDescriptorClass(x.name) for x in self.training.descriptors.all()]
        self.algorithmClass = self.findAlgorithmClass(self.training.algorithm.name)
        self.metricClasses = [self.findMetricClass(x.name) for x in self.validation.metrics.all()]
        self.model = None
        self.onFitCall = onFitCall
        self.molset = self.instance.molset
        self.progress = progress
        self.X = None
        self.y = None
        self.progressStages = []
        self.currentProgress = 0
        self.errors = []

    def recordProgress(self):
        if self.progress and self.currentProgress <= len(self.progressStages):
            self.progress.set_progress(
                self.currentProgress
                , len(self.progressStages)
                , description=self.progressStages[self.currentProgress]
            )
            print(self.progressStages[self.currentProgress])
        self.currentProgress += 1

    def calculateDescriptors(self):
        if self.X.empty:
            self.X = DataFrame() # TODO: implement
            self.recordProgress()

    def saveActivities(self):
        if not self.y:
            activity_set = self.instance.molset.activities.get()
            compounds, activities = activity_set.cleanForModelling()
            self.y = Series(activities) # TODO: implement (take mode into account)
            self.X = DataFrame({"SMILES" : [x.canonicalSMILES for x in compounds]})
            self.recordProgress()

    def fit(self, callback=None) -> models.QSARModel:
        # TODO: check if length of X and y are the same
        self.model = self.algorithmClass(self.training, callback if callback else self.onFitCall)
        self.recordProgress()
        self.model.fit(self.X, self.y)
        return self.instance

    def fitValidate(self) -> models.QSARModel:
        ret = self.fit()
        self.recordProgress()
        self.validate(self.model, self.X, self.y)
        return ret

    def validate(self, model, X : DataFrame, y_truth : Series, perfClass=models.ModelPerformance, **kwargs):
        predictions = model.predict(X)
        for metric_class in self.metricClasses:
            perfClass.objects.create(
                metric=models.ModelPerformanceMetric.objects.get(name=metric_class.name),
                value=metric_class(y_truth, predictions),
                model=self.instance
                **kwargs
            )