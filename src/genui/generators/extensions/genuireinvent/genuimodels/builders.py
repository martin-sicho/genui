# genui/generators/extensions/genuireinvent/genuimodels/builders.py
from genui.models.genuimodels import bases
from genui.models.models import Model
from ..models import ReinventNet

from abc import ABC, abstractmethod

class ReinventBuilder(bases.ProgressMixIn, bases.ModelBuilder, ABC):

    @property
    def corePackage(self):
        from .. import genuimodels
        return genuimodels

    def getY(self):
        return None

    @abstractmethod
    def sample(self, n_samples, from_inputs=None):
        pass

class ReinventNetBuilder(bases.ProgressMixIn, bases.ModelBuilder):
    """
    Keep builder flow identical to DrugEx:
      - progress stages around corpus creation
      - getX returns (train, valid) placeholders (REINVENT uses same file for both)
    """

    def __init__(self, instance: ReinventNet, initial: ReinventNet = None, progress=None):
        super().__init__(instance, progress, None)
        # super().__init__(instance, progress, getattr(instance, "validationStrategy", None))
        self.initial = initial
        self.progressStages.append("Creating Corpus...")
        self.progressStages.append("Corpus Done.")

    def getX(self, update=True):
        # Stage 1: “Creating Corpus…”
        self.recordProgress()

        if update:
            # prepareData() writes the cleaned corpus to disk and syncs preview into AUX ModelFile
            train_mf, valid_mf = self.instance.prepareData()
        else:
            # reuse already-prepared artifacts (both train/valid point to the same corpus)
            train_mf = self.instance.corpusFileTrain
            valid_mf = self.instance.corpusValidFile

        # Stage 2: “Corpus Done.”
        self.recordProgress()

        # Return whatever your Algorithm/Model expects as X.
        # Since our REINVENT "model" is CLI-based, returning the ModelFiles is fine;
        # the algorithm can ignore contents and let the instance supply paths.
        return (train_mf, valid_mf)

    def getY(self):
        return None

    def build(self) -> Model:
        if self.instance.molset and self.validation:
            return super().build()
        raise NotImplementedError(
            "Building Reinvent network without MolSet and validation strategy is not allowed."
        )
