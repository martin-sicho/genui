"""
builders

Created by: Martin Sicho
On: 1/26/20, 6:27 PM
"""
from abc import ABC, abstractmethod

from genui.models.genuimodels import bases
from genui.models.models import Model
from .monitors import DrugExMonitor
from ..models import DrugExNet, DrugExAgent

class DrugExBuilder(bases.ProgressMixIn, bases.ModelBuilder, ABC):

    @property
    def corePackage(self):
        from .. import genuimodels
        return genuimodels

    def getY(self):
        return None

    @abstractmethod
    def sample(self, n_samples, from_inputs=None):
        pass

class DrugExNetBuilder(DrugExBuilder):

    def __init__(self, instance: DrugExNet, initial: DrugExNet =None, progress=None, onFit=None):
        super().__init__(instance, progress, onFit)
        self.initial = initial
        self.onFit = DrugExMonitor(self.instance, lambda epoch : self.recordProgress())

        self.progressStages.append("Creating Corpus...")
        self.progressStages.append("Corpus Done.")

    def getX(self, update=True):
        vocabulary = None
        if self.initial:
            vocabulary = self.initial.corpusTrain.getVoc() + self.initial.corpusTest.getVoc()

        self.recordProgress()
        if update:
            corpusTrain, corpusTest = self.instance.prepareData(vocabulary) # populates corpus data sets
        else:
            corpusTrain, corpusTest = self.instance.corpusTrain, self.instance.corpusTest
        self.recordProgress()

        return corpusTrain, corpusTest

    def build(self) -> Model:
        if self.instance.molset and self.validation:
            return super().build()
        else:
            raise NotImplementedError("Building DrugEx network without molecule set and validation strategy is not allowed.")

    def sample(self, n_samples, from_inputs=None):
        return self.model.sample(n_samples, from_inputs)

class DrugExAgentBuilder(DrugExBuilder):

    def __init__(
            self,
            instance: DrugExAgent,
            progress=None,
            onFit=None
    ):
        super().__init__(instance, progress, onFit)
        self.onFit = DrugExMonitor(self.instance, lambda epoch : self.recordProgress())
        self.exploitNet = self.instance.exploitationNet
        self.exploreNet = self.instance.explorationNet
        self.environ = self.instance.environment

    def getX(self):
        return self.exploreNet.corpusTrain, self.exploreNet.corpusTest

    def sample(self, n_samples, from_inputs=None):
        agent_net_builder = DrugExNetBuilder(self.instance.exploitationNet)
        model = agent_net_builder.model.deserialize(self.instance.modelFile.path)
        return model.sample(n_samples, from_inputs)