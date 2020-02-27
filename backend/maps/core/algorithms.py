"""
algorithms

Created by: Martin Sicho
On: 25-02-20, 15:13
"""
from abc import ABC, abstractmethod

import openTSNE
from pandas import DataFrame, Series

from maps.models import Point
from modelling.core.bases import Algorithm
from modelling.models import ModelParameter

class MapAlgorithm(Algorithm, ABC):

    @abstractmethod
    def getPoints(self) -> [Point]:
        pass

class TSNE(MapAlgorithm):
    name = "TSNE"
    parameters = {
        "n_iter" : ModelParameter.INTEGER,
        "perplexity" : ModelParameter.INTEGER,
        "early_exaggeration" : ModelParameter.INTEGER,
        "early_exaggeration_iter" : ModelParameter.INTEGER,
        "exaggeration" : ModelParameter.INTEGER
    }

    class OpenTSNECallback:

        def __init__(self, builder):
            self.builder = builder

        def __call__(self, iteration, error, embedding):
            print(iteration, error)
            self.builder.recordProgress()
            self.saveModelFile()

        def saveModelFile(self):
            self.builder.saveFile()

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)

        if "n_iter" not in self.params:
            self.params["n_iter"] = 750
        if "early_exaggeration_iter" not in self.params:
            self.params["early_exaggeration_iter"] = 250

        if not self.callback:
            stages = [f"Iteration {25 * (x+1)}" for x in range(int((self.params["n_iter"] + self.params["early_exaggeration_iter"]) / 25))]
            if not stages:
                stages.append("Iteration 1")
            builder.progressStages.extend(stages)
            self.callback = self.OpenTSNECallback(builder)

    @classmethod
    def getModes(cls):
        return [cls.MAP]

    @property
    def model(self):
        return self._model

    def fit(self, X: DataFrame, y: Series):
        tsne = openTSNE.TSNE(
            n_components=2
            , neighbors="exact"
            , negative_gradient_method="bh"
            , callbacks=self.callback
            , callbacks_every_iters=25
            , **self.params
        )
        self._model = tsne.fit(X.values)
        return self.getPoints()

    def predict(self, X: DataFrame):
        return self.model.transform(X)

    def getPoints(self) -> [Point]:
        mols = self.builder.mols
        embedding = self.model

        points = []
        if mols.count() == embedding.shape[0]:
            for idx, mol in enumerate(mols.all()):
                x = embedding[idx, 0]
                y = embedding[idx, 1]
                try:
                    point = Point.objects.get(
                        map=self.builder.instance,
                        molecule=mol,
                    )
                    print("Existing point found. It will be updated:")
                    print(f"x: {point.x} -> {x}")
                    print(f"y: {point.y} -> {y}")
                    point.x = x
                    point.y = y
                    point.save()
                except Point.DoesNotExist:
                    point = Point.objects.create(
                        map=self.builder.instance,
                        molecule=mol,
                        x=x,
                        y=y,
                    )
                points.append(point)
            return points
        else:
            raise Exception("The number of items in the embedding is different from the number of original molecules.")