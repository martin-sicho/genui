import logging
import traceback
import pandas as pd
from pandas import DataFrame, Series
from qsprpred.data import MoleculeTable

from genui.compounds.models import Molecule
from genui.models.genuimodels.bases import PredictionMixIn, ModelBuilder, ProgressMixIn
from genui.qsar.genuimodels.bases import EmbeddingBuilderMixIn
from genui.maps import models

logger = logging.getLogger(__name__)

class MapBuilder(EmbeddingBuilderMixIn, PredictionMixIn, ProgressMixIn, ModelBuilder):
    def __init__(self, instance: models.Map, progress=None, onFit=None):
        super().__init__(instance, progress, onFit)
        self.dataset = None
        self.mols = Molecule.objects.filter(
            providers__in=[x for x in self.instance.molsets.all()]
        )
        self.progressStages.extend(["Calculated embeddings."])

    @property
    def corePackage(self):
        from .. import genuimodels
        return genuimodels

    def getDataset(self) -> DataFrame:
        if self.dataset is None:
            mols = [x for x in self.mols.all()]
            df = pd.DataFrame({"SMILES":[x.canonicalSMILES for x in mols], "mols": mols})
            self.dataset = MoleculeTable(df=df, name=self.instance.name)
            self.dataset.addDescriptors(self.embeddingCalculators)
            self.recordProgress()
        return self.dataset

    def getPoints(self):
        if self.model:
            return self.model.getPoints(self.dataset)

    def validate(self, validation_strategy):
        pass

    def build(self) -> models.Model:
        try:
            self.progressStages.extend(["Saving points...", "Serializing as ChemSpaceJS JSON...", "Done."])
            self.recordProgress()

            self.dataset = self.getDataset()
            _, X = self.model.prepareDataset(self.dataset)
            self.model.fit(X)
            points = self.getPoints()
            if not points:
                raise ValueError("Failed to generate points")
            self.recordProgress()

            try:
                self.instance.saveChemSpaceJSON()
            except Exception as e:
                logger.error(f"Error saving ChemSpaceJS JSON for Map {self.instance.pk}: {str(e)} {traceback.format_exc()}")
                raise

            self.recordProgress()
            return self.instance
        except Exception as e:
            logger.error(f"Error building Map {self.instance.pk}: {str(e)} {traceback.format_exc()}")
            raise