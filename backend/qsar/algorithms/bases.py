"""
algorithms

Created by: Martin Sicho
On: 14-01-20, 10:16
"""
from abc import ABC, abstractmethod
from random import randint

from sklearn.model_selection import KFold, StratifiedKFold

from qsar import models
from pandas import DataFrame, Series


class Algorithm(ABC):
    name = None
    CLASSIFICATION = 'classification'
    REGRESSION = 'regression'
    MODES = [
       (CLASSIFICATION, 'Classification'),
       (REGRESSION, 'Regression'),
    ]

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

    def __init__(
            self,
            instance : models.QSARModel,
            training : models.QSARTrainingStrategy,
            validation : models.ValidationStrategy,
            onFitCall=None
    ):
        self.training = training
        self.validation = validation
        self.descriptors = self.training.descriptors
        self.algorithmClass = self.findAlgorithmClass(self.training.algorithm.name)
        self.metricClasses = [self.findMetricClass(x.name) for x in self.validation.metrics.all()]
        self.model = None
        self.onFitCall = onFitCall
        self.instance = instance
        self.molset = self.instance.molset
        self.X = None
        self.y = None
        self.calculateDescriptors()
        self.saveActivities()

    def calculateDescriptors(self):
        if not self.X:
            self.X = DataFrame() # TODO: implement

    def saveActivities(self):
        if not self.y:
            self.y = Series() # TODO: implement (take mode into account)

    def fit(self, callback=None) -> models.QSARModel:
        # TODO: check if length of X and y are the same
        self.model = self.algorithmClass(self.training, callback if callback else self.onFitCall)
        self.model.fit(self.X, self.y)
        return self.instance

    def fitValidate(self) -> models.QSARModel:
        ret = self.fit()
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



class BasicQSARModelBuilder(QSARModelBuilder):

    def __init__(
            self
            , training: models.QSARTrainingStrategy
            , validation: models.BasicValidationStrategy
            , onFitCall = None
            , onCVFitCall = None
            , onValidFitCall = None
    ):
        super().__init__(training, validation, onFitCall)
        self.onCVFitCall = onCVFitCall
        self.onValidFitCall = onValidFitCall
        random_state = randint(0, 2**32 - 1)
        self.X_valid = self.X.sample(frac=self.validation.validSetSize, random_state=random_state)
        self.X_train = self.X.drop(self.X_valid.index)
        self.y_valid = self.y.sample(frac=self.validation.validSetSize, random_state=random_state)
        self.y_train = self.y.drop(self.y_valid.index)

    def fitValidate(self) -> models.QSARModel:

        is_regression = self.training.mode == models.TrainingStrategy.REGRESSION
        if is_regression:
            folds = KFold(self.validation.cvFolds).split(self.X_train)
        else:
            folds = StratifiedKFold(self.validation.cvFolds).split(self.X_train, self.y_train)
        for i, (trained, validated) in enumerate(folds):
            model = self.algorithmClass(self.training)
            model.fit(self.X_train[trained], self.y_train[trained],)
            self.validate(model, self.X_train[validated], self.y_train[validated], perfClass=models.ModelPerformanceCV, fold=i)
            self.onCVFitCall(self, i)

        model = self.algorithmClass(self.training)
        model.fit(self.X_train, self.y_train)
        self.validate(model, self.X_valid, self.y_valid)
        self.onValidFitCall(self)
        return self.fit()