"""
builders

Created by: Martin Sicho
On: 25-02-20, 15:13
"""
from pandas import DataFrame, Series

from compounds.models import Molecule
from maps.core.algorithms import MapAlgorithm
from modelling.core.bases import PredictionMixIn, ModelBuilder, ProgressMixIn
from qsar.core.bases import DescriptorBuilderMixIn
from maps import models

class MapBuilder(DescriptorBuilderMixIn, PredictionMixIn, ProgressMixIn, ModelBuilder):

    def __init__(self, instance: models.Map, progress=None, onFitCall=None):
        super().__init__(instance, progress, onFitCall)
        self.mols = Molecule.objects.filter(
            providers__in=[x for x in self.instance.molsets.all()]
        )

    def getY(self) -> Series:
        pass

    def getX(self) -> DataFrame:
        self.calculateDescriptors(self.mols.all())
        return self.X

    def getPoints(self):
        if self.model:
            return self.model.getPoints()


