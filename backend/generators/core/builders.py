"""
builders

Created by: Martin Sicho
On: 1/26/20, 6:27 PM
"""
from pandas import DataFrame, Series

from modelling.core import bases
from generators import models

class DrugExBuilder(bases.ModelBuilder):

    def __init__(self, instance: models.DrugExNet):
        super().__init__(instance)

    def getY(self) -> Series:
        pass

    def getX(self) -> DataFrame:
        pass