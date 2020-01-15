"""
algorithms

Created by: Martin Sicho
On: 14-01-20, 10:16
"""
from abc import ABC, abstractmethod

from qsar import models
from pandas import DataFrame, Series


class Algorithm(ABC):
    name = None
    CLASSIFICATION = 'classification'
    REGRESSION = 'regression'

    @staticmethod
    @abstractmethod
    def getDjangoModel() -> models.Algorithm:
        pass

    @staticmethod
    @abstractmethod
    def getParams():
        pass

    @staticmethod
    @abstractmethod
    def getModes():
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
    @abstractmethod
    def getDjangoModel(cls):
        pass

class BaseAlgorithm(Algorithm, ABC):
    name = None

    @staticmethod
    def getModes():
        return [BaseAlgorithm.CLASSIFICATION, BaseAlgorithm.REGRESSION]

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
        for mode in BaseAlgorithm.getModes():
            mode = models.AlgorithmMode.objects.get_or_create(name=mode)[0]
            ret.validModes.add(mode)
        ret.save()
        return ret

class QSARModelBuilder:
    # TODO: load this list dynamically based on introspection of a module where descriptor implementations will reside
    SUPPORTED_DESCRIPTORS = [
        "MORGANFP"
    ]

    @staticmethod
    def findAlgorithmClass(name):
        for class_ in BaseAlgorithm.__subclasses__():
            if hasattr(class_, 'name'):
                if name == class_.name:
                    return class_
            else:
                raise Exception("Unspecified name attribute in algorithm subclass: ", repr(class_))

    @staticmethod
    def findMetricClass(name):
        for class_ in ValidationMetric.__subclasses__():
            if hasattr(class_, 'name'):
                if name == class_.name:
                    return class_
            else:
                raise Exception("Unspecified name attribute in metric subclass: ", repr(class_))

    @classmethod
    def getDescriptorGroupsAsModels(cls):
        ret = []
        for group in cls.SUPPORTED_DESCRIPTORS:
            ret.append(models.DescriptorGroup.objects.get_or_create(name=group))
        return ret

    def __init__(
            self,
            instance : models.QSARModel,
            progress = None,
            onFitCall=None
    ):
        self.training = instance.trainingStrategy
        self.validation = instance.validationStrategy
        self.descriptors = self.training.descriptors
        self.algorithmClass = self.findAlgorithmClass(self.training.algorithm.name)
        self.metricClasses = [self.findMetricClass(x.name) for x in self.validation.metrics.all()]
        self.model = None
        self.onFitCall = onFitCall
        self.instance = instance
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
        if not self.X:
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