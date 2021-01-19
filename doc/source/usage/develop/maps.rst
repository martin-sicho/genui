..  _dev-guide-maps-top:

Chemical Space Maps
===================

The last, but one of the most important features of GenUI, is chemical space visualization.
The `genui.maps` application is responsible for generating 2D representations of compound sets.
In GenUI we call these representations chemical space maps and their calculations is
handled in a similar fashion as QSAR models. The only difference between
generating a chemical space map and a QSAR model is that QSAR models map the
data matrix :command:`X` onto a single array of activity values whereas in the case
of a chemical space map the algorithm is projecting onto a 2D matrix of values.

In this tutorial, we will further extend the :code:`qsarextra` package we have created
before (see :ref:`dev-guide-create-qsar-ext`) so make sure to review that part of the tutorial
before continuing. All we need to do is to define a new algorithm in the
:code:`qsarextra.genuimodels.algorithms` module where we already defined the SVM machine
learning algorithm. However, instead of directly subclassing :py:class:`~genui.models.genuimodels.bases.Algorithm`, we
subclass the `MapAlgorithm`. A simple class implementing all necessary methods could
look like this:

..  code-block:: python

    class MyMap(MapAlgorithm):
        name = 'MyMap'

        def getPoints(self, mols, X) -> [Point]:
            pass

        @property
        def model(self):
            return self._model

        def fit(self, X: DataFrame, y = None):
            pass

        def predict(self, X: DataFrame) -> DataFrame:
            pass

You can see that the difference between :py:class:`~genui.models.genuimodels.bases.Algorithm` and `MapAlgorithm`
is that `MapAlgorithm` defines the :py:meth:`~genui.maps.genuimodels.algorithms.MapAlgorithm.getPoints` method and the return
value from :py:meth:`~genui.maps.genuimodels.algorithms.MapAlgorithm.predict` is not a :code:`Series`, but a :code:`DataFrame`
hinting at the fact that it is a matrix. We will get to the :py:meth:`~genui.maps.genuimodels.algorithms.MapAlgorithm.getPoints`
method shortly, but lets show a real world example first. With the help
of the scikit-learn library, we can implement a simple PCA transformation like so:

..  code-block:: python

    class PCA(MapAlgorithm):
        """
        An example integration of Principle Component Analysis (PCA).
        """

        name = 'PCA'

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

            # initialize scikit-learn components
            from sklearn.decomposition import PCA
            self._model = PCA(n_components=2)
            from sklearn.preprocessing import StandardScaler
            self.scaler = StandardScaler()

        def getPoints(self, mols, X) -> [Point]:
            """
            Converts given molecules to points represented
            as the `Point` Django model class.

            Parameters
            ----------
            mols
                The molecules to convert represented by their respective Django model class (all classes or subclasses of `Molecule`). Should have n_mols items.
            X
                The data matrix of shape [n_mols,n_descriptors].

            Returns
            -------
            points : list
                A list of points in the 2D map represented as `Point` instances.
            """

            transformed_data = self.predict(X)
            points = []
            for idx, mol in enumerate(mols):
                x = transformed_data[idx, 0]
                y = transformed_data[idx, 1]
                point = Point.objects.create(
                    # we can get the map being built from our model builder
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
            """
            Fit the transformation on data X.

            Parameters
            ----------
            X
                The data matrix [n_mols,n_descriptors]. The size and values depend on the chosen set of descriptors.
            y
                Not used in maps.

            """

            self.model.fit(self.scale(X))

        def predict(self, X: DataFrame) -> DataFrame:
            """
            Transform given data.

            Parameters
            ----------
            X
                The data matrix [n_mols,n_descriptors]. The size and values depend on the chosen set of descriptors.

            Returns
            -------
            DataFrame
                2D representation of X [n_mols, 2].
            """

            return self.model.transform(self.scale(X))

        def scale(self, X) -> DataFrame:
            """
            Scale the data matrix (this is required for PCA). Convert
            each variable to a distribution with zero mean and unit variance.

            Parameters
            ----------
            X
                The data matrix [n_mols,n_descriptors]. The size and values depend on the chosen set of descriptors.

            Returns
            -------
            DataFrame
                Scaled matrix of the same shape as X [n_mols,n_descriptors].
            """

            return self.scaler.fit_transform(X)

Without docstrings and comments this class would be quite short. When a new
map is created, a `MapBuilder` instance is initialized inside the Celery task,
just like it happens with `QSARModelBuilder`. The map builder is not much different
from the QSAR model builder:

    1. The `MapBuilder` calls the :py:meth:`~genui.maps.genuimodels.algorithms.MapAlgorithm.fit` method of `MapAlgorithm` and supplies the appropriate data matrix of descriptors calculated for compounds in the chosen compound sets.
    2. The builder calls the :py:meth:`~genui.maps.genuimodels.algorithms.MapAlgorithm.getPoints` method with the same matrix :command:`X` as :py:meth:`~genui.maps.genuimodels.algorithms.MapAlgorithm.fit`. The :py:meth:`~genui.maps.genuimodels.algorithms.MapAlgorithm.getPoints` method is similar to :py:meth:`~genui.maps.genuimodels.algorithms.MapAlgorithm.predict`, but instead of returning a 2D matrix representation of the chemical space map, it creates entries in the database that represent the transformed data as instances of `Point`.
    3. Finally, the builder calls the `saveChemSpaceJSON` method on the `Map` instance sothat it can be visualized with `ChemSpaceJS <https://chemspace.img.cas.cz/>`_. We do not need to care much for this step, but it is good to know about it since the generated JSON file is a fast way to display our map on web pages.
