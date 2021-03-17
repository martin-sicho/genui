from pandas import DataFrame, Series

from genui.maps.genuimodels.algorithms import MapAlgorithm
from genui.models.genuimodels.bases import Algorithm
from genui.maps.models import Point
from genui.models.models import ModelParameter


class MyMap(MapAlgorithm):
    name = 'MyMap'

    def getPoints(self, mols, X) -> [Point]:
        pass

    @property
    def model(self):
        return self._model

    def fit(self, X: DataFrame, y=None):
        pass

    def predict(self, X: DataFrame) -> DataFrame:
        pass

class MDS(MapAlgorithm):

    name = "Multidimensional scaling"

    parameters = {
        "metric": {
            "type": ModelParameter.BOOL
            , "defaultValue": True
        },
        "n_init": {
            "type": ModelParameter.INTEGER
            , "defaultValue": 4
        }
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        from sklearn.manifold import MDS
        self._model = MDS(metric=self.params['metric'],
                         n_init=self.params['n_init'])

    def getPoints(self, mols, X: DataFrame) -> [Point]:
        res = self.predict(X)
        points = []
        for idx, mol in enumerate(mols):
            x = res[idx, 0]
            y = res[idx, 1]
            point = Point.objects.create(
                map=self.builder.instance,
                molecule=mol,
                x=x,
                y=y,
            )
            points.append(point)

        return points

    def fit(self, X: DataFrame, y=None):
        self.model.fit(X)

    def predict(self, X: DataFrame) -> DataFrame:
        return self.model.fit_transform(X);

    @property
    def model(self):
        return self._model;


class PCA(MapAlgorithm):

    name = 'ExamplePCA'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        from sklearn.decomposition import PCA
        self._model = PCA(n_components=2)
        from sklearn.preprocessing import StandardScaler
        self.scaler = StandardScaler()

    def getPoints(self, mols, X) -> [Point]:

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

    @property
    def model(self):
        return self._model

    def fit(self, X: DataFrame, y = None):

        self.model.fit(self.scale(X))

    def predict(self, X: DataFrame) -> DataFrame:

        return self.model.transform(self.scale(X))

    def scale(self, X) -> DataFrame:

        return self.scaler.fit_transform(X)

