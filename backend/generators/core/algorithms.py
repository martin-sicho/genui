"""
algorithms

Created by: Martin Sicho
On: 1/26/20, 5:43 PM
"""
import torch

from drugex.api.corpus import Corpus
from drugex.api.pretrain.generators import BasicGenerator
from drugex.api.pretrain.serialization import GeneratorSerializer, StateProvider
from modelling.core import bases
from modelling.models import ModelParameter, ModelFileFormat


class StateSerializer(StateProvider, GeneratorSerializer):

    def __init__(self, path):
        self.path = path

    def saveGenerator(self, generator):
        state = generator.getState()
        torch.save(state, self.path)
        return self.path

    def getState(self, path=None):
        if not path:
            path = self.path
        return torch.load(path)

class DrugExNetwork(bases.Algorithm):
    name = "DrugExNetwork"
    GENERATOR = 'generator'

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)
        self.corpus = None

    def initSelf(self, X):
        self.corpus = X
        if self.builder.initial:
            self.deserialize(self.builder.initial.modelFile.path)
        else:
            self._model = BasicGenerator(
                monitor=self.callback
                , corpus=self.corpus
                , train_params={
                    "epochs" : self.params['nEpochs']
                    , "monitor_freq" : self.params['monitorFrequency']
                }
            )

    @classmethod
    def getModes(cls):
        return [cls.GENERATOR]

    @classmethod
    def getFileFormats(cls, attach_to=None):
        formats = [ModelFileFormat.objects.get_or_create(
            fileExtension=".torch.pkg",
            description="State of a neural network built with pytorch."
        )[0]]
        if attach_to:
            cls.attachToInstance(attach_to, formats, attach_to.fileFormats)

    @staticmethod
    def getParams():
        names = ['nEpochs', 'monitorFrequency'] # TODO: add batch size
        types = [ModelParameter.INTEGER, ModelParameter.INTEGER]
        assert len(names) == len(types)
        return [
            ModelParameter.objects.get_or_create(name=name, contentType=type_, algorithm=DrugExNetwork.getDjangoModel())[0] for name, type_ in zip(names, types)
        ]

    @property
    def model(self):
        return self._model

    def fit(self, X: Corpus, y=None):
        self.initSelf(X)
        valid_set_size = self.validationInfo.validSetSize
        if valid_set_size:
            self.model.pretrain(validation_size=valid_set_size)
        else:
            self.model.pretrain()

    def predict(self, X: Corpus):
        # TODO: implement generation of compounds
        pass

    def getSerializer(self):
        return lambda filename : self.model.save(StateSerializer(filename))

    def deserialize(self, filename):
        state = StateSerializer(filename)
        if not self.corpus:
            self.corpus = self.builder.getX()
        self._model = BasicGenerator(
            monitor=self.callback
            , initial_state=state
            , corpus=self.corpus
            , train_params={
                "epochs" : self.params['nEpochs']
                , "monitor_freq" : self.params['monitorFrequency']
            }
        )