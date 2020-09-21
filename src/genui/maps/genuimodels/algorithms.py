"""
algorithms

Created by: Martin Sicho
On: 25-02-20, 15:13
"""
from abc import ABC, abstractmethod

import openTSNE
from pandas import DataFrame, Series

from genui.maps.models import Point
from genui.models.genuimodels.bases import Algorithm
from genui.models.models import ModelParameter

class MapAlgorithm(Algorithm, ABC):

    @abstractmethod
    def getPoints(self) -> [Point]:
        pass

class TSNE(MapAlgorithm):
    name = "TSNE"
    parameters = {
        "n_iter" : {
            "type" : ModelParameter.INTEGER
            , "defaultValue" : 500
        },
        "perplexity" : {
            "type" : ModelParameter.INTEGER
            , "defaultValue" : 30
        },
        "early_exaggeration" :{
            "type" : ModelParameter.INTEGER
            , "defaultValue" : 12
        },
        "early_exaggeration_iter" : {
            "type" : ModelParameter.INTEGER
            , "defaultValue" : 250
        },
        "exaggeration" : {
            "type" : ModelParameter.INTEGER
            , "defaultValue" : 0
        }
    }

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)

        if "n_iter" not in self.params:
            self.params["n_iter"] = 750
        if "early_exaggeration_iter" not in self.params:
            self.params["early_exaggeration_iter"] = 250
        if "exaggeration" in self.params and self.params["exaggeration"] == 0:
            del self.params["exaggeration"]

        stages = [f"Iteration {25 * (x+1)}" for x in range(int((self.params["n_iter"] + self.params["early_exaggeration_iter"]) / 25))]
        if not stages:
            stages.append("Iteration 1")
        self.builder.progressStages.extend(stages)

    @classmethod
    def getModes(cls):
        return [cls.MAP]

    @property
    def model(self):
        return self._model

    def fit(self, X: DataFrame, y: Series):
        def callback(iteration, error, embedding):
            print(f'Current Iteration: {iteration} (error: {error})')
            self.builder.recordProgress()

        tsne = openTSNE.TSNE(
            n_components=2
            , neighbors="exact"
            , negative_gradient_method="bh"
            , callbacks=callback
            , callbacks_every_iters=25
            , **self.params
        )
        self._model = tsne.fit(X.values)
        return self.getPoints()

    def predict(self, X: DataFrame):
        return self.model.transform(X)

    def getSerializer(self):
        del self.model.gradient_descent_params['callbacks']
        return super().getSerializer()

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