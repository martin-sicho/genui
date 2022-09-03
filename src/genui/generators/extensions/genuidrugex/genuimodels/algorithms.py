"""
algorithms

Created by: Martin Sicho
On: 1/26/20, 5:43 PM
"""
from abc import ABC
import torch

from genui.models.genuimodels import bases
from genui.models.models import ModelParameter, ModelFileFormat

class DrugExAlgorithm(bases.Algorithm, ABC):

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)
        if 'nEpochs' in self.params:
            for i in range(self.params['nEpochs']):
                self.builder.progressStages.append(f"Epoch {i+1}")
        self.builder.progressStages.append("Model Built.")
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

    def predict(self, X):
        return self.model.sample(X)

    def sample(self, n_samples, from_inputs=None):
        batch_size = min(self.params['batchSize'] if 'batchSize' in self.params else 32, n_samples) # FIXME: move this up from subclasses
        if from_inputs:
            inputs = from_inputs.asDataLoader(batch_size)
        else:
            inputs = self.builder.getX(update=False)[0].asDataLoader(batch_size)
        if batch_size == n_samples:
            for x in inputs:
                inputs = [x]
                break
        else:
            inputs = [x for i,x in enumerate(inputs) if i <= n_samples // batch_size]
        smiles, frags = self.predict(inputs)
        return smiles, frags

    def getSerializer(self):
        return lambda path : torch.save(self.model.getModel(), path)

    def getDeserializer(self):
        def deserializer(path):
            self.model.loadStatesFromFile(path)
            return self.model
        return deserializer

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
        }
    }

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)
        self._model = self.builder.instance.getModel()

    def fit(self, X, y=None):
        self._model = self.builder.instance.getModel()
        if self.builder.initial:
            # load initial states if finetuning
            self.deserialize(self.builder.initial.modelFile.path)
        self._model.fit(
            X[0].asDataLoader(self.params['batchSize']),
            X[1].asDataLoader(self.params['batchSize']),
            epochs=self.params['nEpochs'],
            monitor=self.callback
        )

class DrugExAgent(DrugExAlgorithm):
    name = "DrugExAgent"
    parameters = {
        'nEpochs': {
            "type" : ModelParameter.INTEGER,
            "defaultValue" : 60
        },
        'batchSize' : {
            "type" : ModelParameter.INTEGER,
            "defaultValue" : 512
        },
        'epsilon' : {
            "type" : ModelParameter.FLOAT,
            "defaultValue" : 0.01
        },
        'beta' : {
            "type" : ModelParameter.FLOAT,
            "defaultValue" : 0.1
        }
    }

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)
        self.environ = self.instance.environment.getInstance()
        self.exploitNet = self.instance.exploitationNet.getModel()
        self.exploreNet = self.instance.explorationNet.getModel()
        self._model = self.instance.trainingStrategy.getExplorerInstance(self.exploitNet, self.environ, self.exploreNet, self.params['epsilon'], self.params['beta'], self.params['batchSize'])

    def fit(self, X=None, y=None):
        self.model.fit(
            X[0].asDataLoader(self.params['batchSize']),
            X[1].asDataLoader(self.params['batchSize']),
            monitor=self.callback,
            epochs=self.params['nEpochs']
        )
        self.builder.recordProgress()