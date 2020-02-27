"""
builders

Created by: Martin Sicho
On: 25-02-20, 15:13
"""
from pandas import DataFrame, Series

from compounds.models import Molecule
from modelling.core.bases import PredictionMixIn, ModelBuilder, ProgressMixIn
from qsar.core.bases import DescriptorBuilderMixIn
from maps import models

class MapBuilder(DescriptorBuilderMixIn, PredictionMixIn, ProgressMixIn, ModelBuilder):

    def __init__(self, instance: models.Map, progress=None, onFit=None):
        super().__init__(instance, progress, onFit)
        self.mols = Molecule.objects.filter(
            providers__in=[x for x in self.instance.molsets.all()]
        )
        self.progressStages.extend(["Initialized.", "Calculated descriptors."])
        self.recordProgress()

    def getY(self) -> Series:
        pass

    def getX(self) -> DataFrame:
        self.calculateDescriptors(self.mols.all())
        self.recordProgress()
        return self.X

    def getPoints(self):
        if self.model:
            return self.model.getPoints()


