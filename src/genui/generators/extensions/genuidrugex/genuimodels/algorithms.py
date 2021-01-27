"""
algorithms

Created by: Martin Sicho
On: 1/26/20, 5:43 PM
"""
from abc import ABC

import torch
import pandas as pd

from drugex.api.agent.agents import DrugExAgent as DrugExAgentTrainer
from drugex.api.agent.policy import PG
from drugex.api.corpus import Corpus, DataProvidingCorpus
from drugex.api.environ.models import EnvironProvider
from drugex.api.pretrain.generators import BasicGenerator, Generator
from drugex.api.pretrain.serialization import GeneratorSerializer, StateProvider, GeneratorDeserializer
from genui.generators.extensions.genuidrugex.models import DrugExNet
from genui.models.genuimodels import bases
from genui.models.models import ModelParameter, ModelFileFormat
from genui.qsar.genuimodels.bases import QSARModelBuilder
from genui.qsar.models import QSARModel

class StateSerializer(StateProvider, GeneratorDeserializer, GeneratorSerializer):

    def getGenerator(self):
        if not self.corpus:
            raise Exception("Corpus must be specified for deserialization.")
        return BasicGenerator(
            monitor=self.monitor
            , initial_state=self
            , corpus=self.corpus
            , train_params=self.train_params
        )

    def __init__(self, state_path, corpus=None, monitor=None, train_params=None):
        self.statePath = state_path
        self.corpus = corpus
        self.monitor = monitor
        self.train_params = train_params

    def saveGenerator(self, generator):
        state = generator.getState()
        torch.save(state, self.statePath)
        return self.statePath

    def getState(self, path=None):
        if not path:
            path = self.statePath
        if torch.cuda.is_available():
            return torch.load(path, map_location=torch.device('cuda', torch.cuda.current_device()))
        else:
            return torch.load(path, map_location=torch.device('cpu'))

    @staticmethod
    def getFromModel(model : DrugExNet, monitor=None, train_params=None) -> Generator:
        state = StateSerializer(
            model.modelFile.path,
            model.corpus,
            monitor,
            train_params
        )
        return state.getGenerator()

class EnvironWrapper(EnvironProvider):

    def __init__(self, model : QSARModel):
        self.builder =  QSARModelBuilder(model)

    def predictSMILES(self, smiles):
        self.builder.calculateDescriptors(smiles)
        return self.builder.predict()

class DrugExAlgorithm(bases.Algorithm, ABC):

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)
        if 'nEpochs' in self.params:
            for i in range(self.params['nEpochs']):
                self.builder.progressStages.append(f"Epoch {i+1}")
        self.builder.progressStages.append("Model Built.")
        self.corpus = self.builder.getX()
        self.train_params = dict()

    @classmethod
    def getFileFormats(cls, attach_to=None):
        formats = [
            ModelFileFormat.objects.get_or_create(
            fileExtension=".torch.pkg",
            description="State of a neural network built with pytorch."
        )[0],
        ModelFileFormat.objects.get_or_create(
            fileExtension=".pkg",
            description="State of a neural network built with pytorch."
        )[0]
        ]
        if attach_to:
            cls.attachToInstance(attach_to, formats, attach_to.fileFormats)

    @classmethod
    def getModes(cls):
        return [cls.GENERATOR]

    @property
    def model(self):
        return self._model

    def predict(self, X) -> pd.Series:
        return self.model.sample(X)

    def sample(self, n_samples):
        return self.predict(n_samples)

    def getSerializer(self):
        return lambda filename : StateSerializer(filename).saveGenerator(self.model)

    def getDeserializer(self):
        return lambda filename : StateSerializer(
            filename,
            self.corpus,
            self.callback,
            train_params=self.train_params
        ).getGenerator()

class DrugExNetwork(DrugExAlgorithm):
    name = "DrugExNetwork"
    parameters = {
        'nEpochs': {
            "type" : ModelParameter.INTEGER,
            "defaultValue" : 60
        },
        'batchSize' : {
            "type" : ModelParameter.INTEGER,
            "defaultValue" : 512
        },
        'monitorFrequency' : {
            "type" : ModelParameter.INTEGER,
            "defaultValue" : 10
        },
    }

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)
        self.train_params = {
            "epochs" : self.params['nEpochs'] if 'nEpochs' in self.params else 100
            , "monitor_freq" : self.params['monitorFrequency'] if 'monitorFrequency' in self.params else 100
        }
        self.loaders_params = {"batch_size" : self.params['batchSize'] if 'batchSize' in self.params else 512}

    def initSelf(self, X):
        self.corpus = X
        if self.builder.initial:
            self.deserialize(self.builder.initial.modelFile.path)
        else:
            self._model = BasicGenerator(
                monitor=self.callback
                , corpus=self.corpus
                , train_params=self.train_params
            )

    def fit(self, X: Corpus, y=None):
        if not isinstance(X, DataProvidingCorpus):
            raise NotImplementedError(f"You need an instance of {DataProvidingCorpus.__name__} to fit a DrugEx network.")
        self.initSelf(X)
        valid_set_size = self.validationInfo.validSetSize if self.validationInfo else 0
        tlp = self.loaders_params
        vlp = self.loaders_params
        if valid_set_size:
            self.model.pretrain(validation_size=valid_set_size, train_loader_params=tlp, valid_loader_params=vlp)
        else:
            self.model.pretrain(train_loader_params=tlp, valid_loader_params=vlp)
        self.builder.recordProgress()

class DrugExAgent(DrugExAlgorithm):
    name = "DrugExAgent"
    parameters = {
        'nEpochs': {
            "type" : ModelParameter.INTEGER,
            "defaultValue" : 60
        },
        'pg_batch_size' : {
            "type" : ModelParameter.INTEGER,
            "defaultValue" : 512
        },
        'pg_mc' : {
            "type" : ModelParameter.INTEGER,
            "defaultValue" : 5
        },
        'pg_epsilon' : {
            "type" : ModelParameter.FLOAT,
            "defaultValue" : 0.01
        },
        'pg_beta' : {
            "type" : ModelParameter.FLOAT,
            "defaultValue" : 0.1
        }
    }

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)
        self.environ = self.instance.environment
        self.exploitNet = self.instance.exploitationNet
        self.exploreNet = self.instance.explorationNet
        self._model = DrugExAgentTrainer(
            self.callback # our monitor
            , self.wrapQSARModelToEnviron(self.environ)
            , StateSerializer.getFromModel(self.exploitNet)
            , PG(
                batch_size=self.params['pg_batch_size']
                , mc=self.params['pg_mc']
                , epsilon=self.params['pg_epsilon']
                , beta=self.params['pg_beta']
            )
            , StateSerializer.getFromModel(self.exploreNet)
            , {
                "n_epochs" : self.params['nEpochs']
            }
        )

    @staticmethod
    def wrapQSARModelToEnviron(model : QSARModel):
        return EnvironWrapper(model)

    def fit(self, X=None, y=None):
        self.model.train()
        self.builder.recordProgress()