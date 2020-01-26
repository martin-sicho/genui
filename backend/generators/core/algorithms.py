"""
algorithms

Created by: Martin Sicho
On: 1/26/20, 5:43 PM
"""
from pandas import DataFrame, Series

from modelling.core import bases

class DrugExNetwork(bases.Algorithm):

    def __init__(self, builder):
        super().__init__(builder)

    @staticmethod
    def getParams():
        pass

    @property
    def model(self):
        pass

    def fit(self, X: DataFrame, y: Series):
        pass

    def predict(self, X: DataFrame):
        pass