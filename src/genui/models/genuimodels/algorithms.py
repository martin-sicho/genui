"""
algorithms

Created by: Martin Sicho
On: 24-01-20, 15:03
"""
from pandas import DataFrame, Series
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier

from . import bases
from genui.models.models import ModelParameter


class RandomForest(bases.Algorithm):
    name = "RandomForest"
    parameters = {
        "n_estimators" : {
            "type" : ModelParameter.INTEGER,
            "defaultValue" : 100
        }
    }

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)
        self.alg = RandomForestRegressor if self.mode.name == self.REGRESSION else RandomForestClassifier

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