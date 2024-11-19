import logging
from pandas import DataFrame, Series
from genui.compounds.models import Molecule
from genui.models.genuimodels.bases import PredictionMixIn, ModelBuilder, ProgressMixIn
from genui.qsar.genuimodels.bases import DescriptorBuilderMixIn
from genui.maps import models

logger = logging.getLogger(__name__)

class MapBuilder(DescriptorBuilderMixIn, PredictionMixIn, ProgressMixIn, ModelBuilder):
    def __init__(self, instance: models.Map, progress=None, onFit=None):
        super().__init__(instance, progress, onFit)
        self.mols = Molecule.objects.filter(
            providers__in=[x for x in self.instance.molsets.all()]
        )
        self.progressStages.extend(["Calculated descriptors."])

    @property
    def corePackage(self):
        from .. import genuimodels
        return genuimodels

    def getY(self) -> Series:
        pass

    def getX(self) -> DataFrame:
        if self.X is None:
            self.X = self.calculateDescriptors(self.mols.all())
            self.recordProgress()
        return self.X

    def getPoints(self):
        if self.model:
            # TODO: check that number of mols and rows of X are the same
            return self.model.getPoints(self.mols.all(), self.getX())

    def build(self) -> models.Model:
        try:
            super().build()
            self.progressStages.extend(["Saving points...", "Serializing as ChemSpaceJS JSON...", "Done."])
            self.recordProgress()

            points = self.getPoints()
            if not points:
                raise ValueError("Failed to generate points")
            self.recordProgress()

            try:
                self.instance.saveChemSpaceJSON()
            except Exception as e:
                logger.error(f"Error saving ChemSpaceJS JSON for Map {self.instance.pk}: {str(e)}")
                raise

            self.recordProgress()
            return self.instance
        except Exception as e:
            logger.error(f"Error building Map {self.instance.pk}: {str(e)}")
            raise