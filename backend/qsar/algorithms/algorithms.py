"""
algorithms

Created by: Martin Sicho
On: 14-01-20, 11:05
"""
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

from qsar import models
from . import bases
from pandas import DataFrame, Series


class RandomForest(bases.BaseAlgorithm):
    name = "RandomForest"

    def __init__(self, training_info: models.TrainingStrategy, callback=None):
        super().__init__(training_info, callback)
        self.alg = RandomForestRegressor if self.mode == models.TrainingStrategy.REGRESSION else RandomForestClassifier
        self._model = None

    @staticmethod
    def getParams():
        names = ['n_estimators',]
        types = [models.ModelParameter.INTEGER]
        return [
            models.ModelParameter.objects.get_or_create(name=name, contentType=type_, algorithm=RandomForest.getDjangoModel())[0] for name, type_ in zip(names, types)
        ]

    @staticmethod
    def getDjangoModel() -> models.Algorithm:
        ret = models.Algorithm.objects.get_or_create(
            name=RandomForest.name
        )[0]
        file_format = models.ModelFileFormat.objects.get_or_create(
            fileExtension=".joblib.gz",
            description="A compressed joblib file."
        )[0]
        ret.fileFormats.add(file_format)
        ret.save()
        return ret

    @property
    def model(self):
        return self._model

    def fit(self, X : DataFrame, y : Series):
        self._model = self.alg(**self.params)
        self._model.fit(X, y)
        if self.callback:
            self.callback(self)

    def predict(self, X : DataFrame):
        is_regression = self.trainingInfo.mode == models.TrainingStrategy.REGRESSION
        if self.model:
            if is_regression:
                return self.model.predict(X)
            else:
                return self.model.predict_proba(X)
        else:
            raise Exception("You have to fit the model first.")