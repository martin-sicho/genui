"""
builders

Created by: Martin Sicho
On: 25-02-20, 15:13
"""
from pandas import DataFrame, Series

from genui.compounds.models import Molecule
from genui.modelling.genuimodels.bases import PredictionMixIn, ModelBuilder, ProgressMixIn
from genui.qsar.genuimodels.bases import DescriptorBuilderMixIn
from genui.maps import models

class MapBuilder(DescriptorBuilderMixIn, PredictionMixIn, ProgressMixIn, ModelBuilder):

    def __init__(self, instance: models.Map, progress=None, onFit=None):
        super().__init__(instance, progress, onFit)
        self.mols = Molecule.objects.filter(
            providers__in=[x for x in self.instance.molsets.all()]
        )
        self.progressStages.extend(["Calculated descriptors."])

    @property
    def coreModule(self):
        from .. import genuimodels
        return genuimodels

    def getY(self) -> Series:
        pass

    def getX(self) -> DataFrame:
        self.calculateDescriptors(self.mols.all())
        self.recordProgress()
        return self.X

    def getPoints(self):
        if self.model:
            return self.model.getPoints()


