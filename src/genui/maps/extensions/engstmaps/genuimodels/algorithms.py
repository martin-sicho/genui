from pandas import DataFrame
from qsprpred.data import MoleculeTable
import numpy as np

from genui.maps.genuimodels.algorithms import MapAlgorithm
from genui.maps.models import Point
from genui.models.models import ModelParameter

import umap

class EngstMaps(MapAlgorithm):
    def getPoints(self, dataset: MoleculeTable) -> [Point]:
        mols, X = self.prepareDataset(dataset)
        transformed_data = self.predict(X)
        points = []
        for idx, mol in enumerate(mols):
            x = transformed_data[idx, 0]
            y = transformed_data[idx, 1]
            point = Point.objects.create(
                map=self.builder.instance,
                molecule=mol,
                x=x,
                y=y,
            )
            points.append(point)

        return points


class Isomap(EngstMaps):
    name = "Isomap"
    parameters = {
        "n_neighbors": {
            "type": ModelParameter.INTEGER
            , "defaultValue": 5
        }
    }

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)

        if "n_neighbors" not in self.params:
            self.params["n_neighbors"] = 5
        from sklearn.manifold import Isomap
        self._model = Isomap(n_neighbors=self.params['n_neighbors'])




    def fit(self, X: DataFrame, y=None):
        self._model = self._model.fit(X)

    def predict(self, X: DataFrame) -> DataFrame:
        return self.model.transform(X)

    @property
    def model(self):
        return self._model

class SVD(EngstMaps):
    name = "SVD"

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)

        from sklearn.decomposition import TruncatedSVD
        self._model = TruncatedSVD()

    def fit(self, X: DataFrame, y=None):
        self._model = self._model.fit(X)

    def predict(self, X: DataFrame) -> DataFrame:
        return self.model.transform(X)

    @property
    def model(self):
        return self._model

class UMAP(EngstMaps):
    name = "UMAP"
    parameters = {
        "n_neighbors": {
            "type": ModelParameter.INTEGER
            , "defaultValue": 15
        },
        "min_dist": {
            "type": ModelParameter.FLOAT
            , "defaultValue": 0.1
        },
        "metric": {
            "type": ModelParameter.STRING
            , "defaultValue": 'euclidean'
        }
    }

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)

        if "n_neighbors" not in self.params:
            self.params["n_neighbors"] = 15
        if "min_dist" not in self.params:
            self.params["min_dist"] = 0.1
        if "metric" not in self.params:
            self.params["metric"] = 'euclidean'

        self._model = umap.UMAP(n_neighbors=self.params["n_neighbors"],
                                min_dist=self.params["min_dist"],
                                metric=self.params["metric"])


    def fit(self, X: DataFrame, y=None):
        self._model = self._model.fit(X)

    def predict(self, X: DataFrame) -> DataFrame:
        return self.model.transform(X)

    @property
    def model(self):
        return self._model

class PCA(EngstMaps):

    name = 'PCA'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        from sklearn.decomposition import PCA
        self._model = PCA(n_components=2)
        from sklearn.preprocessing import StandardScaler
        self.scaler = StandardScaler()


    @property
    def model(self):
        return self._model

    def fit(self, X: DataFrame, y = None):

        self._model = self.model.fit(self.scale(X))

    def predict(self, X: DataFrame) -> DataFrame:

        return self.model.transform(self.scale(X))

    def scale(self, X) -> DataFrame:
        X = np.array(X)
        if not X.flags["WRITEABLE"]:
            X = np.copy(X)
        return self.scaler.fit_transform(X)

class MDS(EngstMaps):

    name = "MDS"

    parameters = {
        "n_init": {
            "type": ModelParameter.INTEGER
            , "defaultValue": 4
        }
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if "n_init" not in self.params:
            self.params["n_init"] = 4

        from sklearn.manifold import MDS
        self._model = MDS(n_init=self.params['n_init'])

    def fit(self, X: DataFrame, y=None):
        self._model = self.model.fit(X)

    def predict(self, X: DataFrame) -> DataFrame:
        return self.model.fit_transform(X);

    @property
    def model(self):
        return self._model

