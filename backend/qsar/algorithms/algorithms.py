"""
algorithms

Created by: Martin Sicho
On: 14-01-20, 11:05
"""
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

from qsar import models
from . import bases
from pandas import DataFrame, Series


class RandomForest(bases.Algorithm):
    name = "RandomForest"

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)
        self.alg = RandomForestRegressor if self.mode.name == self.REGRESSION else RandomForestClassifier

    @staticmethod
    def getParams():
        names = ['n_estimators',]
        types = [models.ModelParameter.INTEGER]
        return [
            models.ModelParameter.objects.get_or_create(name=name, contentType=type_, algorithm=RandomForest.getDjangoModel())[0] for name, type_ in zip(names, types)
        ]

    @property
    def model(self):
        return self._model

    def fit(self, X : DataFrame, y : Series):
        self._model = self.alg(**self.params)
        self._model.fit(X, y)
        if self.callback:
            self.callback(self)

    def predict(self, X : DataFrame):
        is_regression = self.trainingInfo.mode.name == self.REGRESSION
        if self.model:
            if is_regression:
                return self.model.predict(X)
            else:
                return self.model.predict_proba(X)
        else:
            raise Exception("You have to fit the model first.")