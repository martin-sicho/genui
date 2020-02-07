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
from drugex.api.corpus import Corpus, CorpusCSV
from drugex.api.environ.models import EnvironProvider
from drugex.api.pretrain.generators import BasicGenerator, Generator
from drugex.api.pretrain.serialization import GeneratorSerializer, StateProvider, GeneratorDeserializer
from generators import models
from modelling.core import bases
from modelling.models import ModelParameter, ModelFileFormat
from qsar.core.bases import QSARModelBuilder
from qsar.models import QSARModel


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
        return torch.load(path)

    @staticmethod
    def getFromModel(model : models.DrugExNet, monitor=None, train_params=None) -> Generator:
        state = StateSerializer(
            model.modelFile,
            CorpusCSV.fromFiles(model.corpus.corpusFile.path, model.corpus.vocFile.path),
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
    GENERATOR = 'generator'

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)
        for i in range(self.params['nEpochs']):
            self.builder.progressStages.append(f"Epoch {i+1}")
        self.corpus = None

    @classmethod
    def getFileFormats(cls, attach_to=None):
        formats = [ModelFileFormat.objects.get_or_create(
            fileExtension=".torch.pkg",
            description="State of a neural network built with pytorch."
        )[0]]
        if attach_to:
            cls.attachToInstance(attach_to, formats, attach_to.fileFormats)

    @classmethod
    def getModes(cls):
        return [cls.GENERATOR]

    @property
    def model(self):
        return self._model

    def predict(self, X) -> pd.Series:
        # TODO: implement generation of compounds
        n_samples = X
        print(n_samples)
        raise NotImplementedError("This is not working yet...")

    def sample(self, n_samples):
        return self.predict(n_samples)

    def getSerializer(self):
        return lambda filename : StateSerializer(filename).saveGenerator(self.model)

    # def getSerializer(self):
    #     return lambda filename : self.model.save(StateSerializer(filename))

    def getDeserializer(self):
        return lambda filename : StateSerializer(filename, self.corpus, self.callback).getGenerator()

class DrugExNetwork(DrugExAlgorithm):
    name = "DrugExNetwork"
    parameters = {
        'nEpochs': ModelParameter.INTEGER,
        'monitorFrequency' : ModelParameter.INTEGER,
    }

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)
        self.train_params = {
            "epochs" : self.params['nEpochs']
            , "monitor_freq" : self.params['monitorFrequency']
        }

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
        self.initSelf(X)
        valid_set_size = self.validationInfo.validSetSize
        if valid_set_size:
            self.model.pretrain(validation_size=valid_set_size)
        else:
            self.model.pretrain()

    def getDeserializer(self):
        if not self.corpus:
            self.corpus = self.builder.getX()
        return lambda filename : StateSerializer(
            filename,
            self.corpus,
            self.callback,
            train_params=self.train_params
        ).getGenerator()

class DrugExAgent(DrugExAlgorithm):
    name = "DrugExAgent"
    parameters = {
        'nEpochs': ModelParameter.INTEGER,
        'pg_batch_size' : ModelParameter.INTEGER,
        'pg_mc' : ModelParameter.INTEGER,
        'pg_epsilon' : ModelParameter.FLOAT,
        'pg_beta' : ModelParameter.FLOAT
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