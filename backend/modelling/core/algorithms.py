"""
algorithms

Created by: Martin Sicho
On: 24-01-20, 15:03
"""
from pandas import DataFrame, Series
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier

from . import bases
import modelling.models


class RandomForest(bases.Algorithm):
    name = "RandomForest"

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)
        self.alg = RandomForestRegressor if self.mode.name == self.REGRESSION else RandomForestClassifier

    @staticmethod
    def getParams():
        names = ['n_estimators',]
        types = [modelling.models.ModelParameter.INTEGER]
        return [
            modelling.models.ModelParameter.objects.get_or_create(name=name, contentType=type_, algorithm=RandomForest.getDjangoModel())[0] for name, type_ in zip(names, types)
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
                return self.model.predict_proba(X)[:,1]
        else:
            raise Exception("You have to fit the model first.") # TODO: add custom exception