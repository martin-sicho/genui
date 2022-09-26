
..  _dev-guide-create-qsar-ext:

QSAR
====

Adding a Model Algorithm
------------------------

In this tutorial we will do a direct follow up of the example we showed in
:ref:`dev-guide-create-extension`. We will add a new machine learning algorithm
for QSAR modelling to the :code:`qsarextra` extension that
we created previously. You will see that we will actually not need to
explicitly define any REST API endpoints or handle web requests in order to
make our implementation visible in the REST API. The GenUI framework
already defines suitable endpoints and will simply call your code when
required. Therefore, we can mainly focus on the implementation of our algorithm.

In order to add new QSAR model implementations, all we have to do is create a new subpackage in our :code:`qsarextra` application. This package
will be called :code:`genuimodels` and its presence tells the :code:`genuisetup`
command that it should look for special modules and classes in this package. For new machine learning algorithms, we need to create a module called
:code:`algorithms.py`. The :code:`genuisetup` command will be looking for a module named
like this and search for implementations of the `genui.models.genuimodels.bases.Algorithm`
abstract class. A minimal :code:`algorithms.py` would look something like this:

..  code-block:: python

    """
    algorithms.py in src/genui/qsar/extensions/qsarextra/genuimodels/

    """

    from pandas import DataFrame, Series
    from genui.models.genuimodels.bases import Algorithm

    class MyAlgorithm(Algorithm):
        name = "MyAlgorithmName"

        @property
        def model(self):
            pass

        def fit(self, X: DataFrame, y: Series):
            pass

        def predict(self, X: DataFrame) -> Series:
            pass

Implement these three methods and you are done. No more work needed.
The algorithm should now show up among the others in the REST API
(URL: :code:`/api/qsar/algorithms/`)
after you run the :code:`genuisetup` command.

What happens after you run :code:`genuisetup`
is that the new algorithm will be registered and an entry will be created in the
database representing this class. If the user then selects this as
an option while defining a QSAR model with the API (or in the GUI), the GenUI framework knows it needs
to use this class to construct the model. It prepares all the data (depending
on what descriptors were chosen by the user) and after initialization the
`Algorithm.fit` method is called. Similarly, for predictions the
`Algorithm.predict` method is used.

Lets see how a real world example would look like. We will implement
a new algorithm based on the implementation of Support Vector Machines (SVMs)
in scikit-learn. We could include SVMs in :code:`qsarextra` by defining
the following class:

..  code-block:: python


    """
    algorithms.py in src/genui/qsar/extensions/qsarextra/genuimodels/

    """

    from pandas import DataFrame, Series
    from sklearn.svm import SVR, SVC

    from genui.models.genuimodels.bases import Algorithm
    from genui.models.models import ModelParameter


    class SVM(Algorithm):
        name = "SVM"
        parameters = {
            "C" : {
                "type" : ModelParameter.FLOAT,
                "defaultValue" : 1.0
            },
            "kernel" : {
                "type" : ModelParameter.STRING,
                "defaultValue" : 'rbf'
            }
        }

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs) # call base class constructor
            self.alg = SVR if self.mode.name == self.REGRESSION else SVC # based on prediction mode, get the correct scikit-learn class

        @property
        def model(self):
            """
            You define this property so that it returns the final fitted model.
            It can be any object so it is ok if we just return the SVC/SVR instance
            directly.

            This object is used mainly for serialization to disk and you can
            implement methods that do the job. GenUI uses *joblib* by default,
            which can handle scikit-learn instances just fine so there
            is no need to customize anything here.

            Returns
            -------
            object
                An instance representing the fitted model.
            """

            return self._model # None by default

        def fit(self, X: DataFrame, y: Series):
            """
            This method takes the data matrix and fits the model.
            The input will be a `DataFrame` and `Series`.
            Data will usually be raw without any transformations
            or normalizations applied so you might want to do them
            here as well.

            Parameters
            ----------
            X : DataFrame
                The data matrix to fit by the model. Samples as rows, variables as columns.
            y : Series
                The ground truth value for each sample. Should be the same length as rows of X.
            """

            # we also want probabilities for classification (see the 'predict' method)
            # so we add the 'probability' parameter when needed
            self._model = self.alg(probability=True, **self.params) if self.alg.__name__  == SVC.__name__ else self.alg(**self.params)

            self._model.fit(X, y)
            if self.callback:
                self.callback(self)

        def predict(self, X: DataFrame) -> Series:
            """
            A method used for predictions. You get
            a matrix of samples (you should again transform
            and normalize and needed) and it is expected
            your model returns the predictions as a `Series`.

            Parameters
            ----------
            X : DataFrame
                The samples.

            Returns
            -------
            predictions : Series
                The predictions.

            """

            is_regression = self.mode.name == self.REGRESSION
            if self.model:
                if is_regression:
                    return self.model.predict(X)
                else:
                    return self.model.predict_proba(X)[:,1]
            else:
                raise Exception("You have to fit the model first.")

For more information on other useful attributes and methods,
see the `genui.models.genuimodels.bases.Algorithm` reference.

Writing Tests
~~~~~~~~~~~~~

It is always good practice to validate newly implemented features with unit tests.
The GenUI framework defines a few classes that make writing tests easier. In order
to test our SVM models, we could define the following test case in the
:code:`qsarextra.tests` module:

..  code-block:: python

    """
    tests.py in src/genui/qsar/extensions/qsarextra/

    """

    from rest_framework.test import APITestCase

    from genui.models.models import AlgorithmMode, Algorithm
    from genui.qsar.tests import QSARModelInit


    class QSARExtraTestCase(QSARModelInit, APITestCase):

        def test_my_SVC(self):
            self.createTestQSARModel(
                mode = AlgorithmMode.objects.get(name="classification"),
                algorithm = Algorithm.objects.get(name="SVM"),
                parameters={
                    "C" : 1.5,
                    "kernel" : 'poly'
                }
            )

        def test_my_SVR(self):
            self.createTestQSARModel(
                mode = AlgorithmMode.objects.get(name="regression"),
                algorithm = Algorithm.objects.get(name="SVM"),
                parameters={
                    "C" : 1.5,
                    "kernel" : 'poly'
                }
            )

The `createTestQSARModel` method of `QSARModelInit` defines a basic unit test
to train a given QSAR model using the REST API. It automatically sets up a project and imports
some test compounds and bioactivites from the ChEMBL database for training.
The resulting model is returned from the method as the appropriate Django model.

..  note:: You can run all tests for GenUI with :code:`python manage.py test`.
    However, you will need to set the settings module to `genui.settings.test`.
    This is the same as the `genui.settings.debug` configuration, but all Celery tasks will be ran
    synchronously in a single thread and created media files are saved into a separate directory while executing tests as well.

Adding New Molecular Descriptors
--------------------------------

In QSAR modelling, an important decision is the choice of molecular descriptors
so you will likely want to implement calculation of your own. Doing so
is easy and it is again done through the definition of a special class.
This time we will need to implement the :code:`DescriptorCalculator.__call__` method
of the `DescriptorCalculator` abstract class defined in the `genui.qsar` package.

Lets say we would like to have the :code:`qsarextra` extension provide
a new set of chemical descriptors. We have to create a new module under
:code:`genui.qsar.extensions.qsarextra.genuimodels`,
but this time we will name it :code:`descriptors.py`.
In this file, we can define the descriptor calculators.
For example, we could include the 2D descriptors provided
by the RDKit library like so:

..  code-block::  python

    """
    descriptors.py in src/genui/qsar/extensions/qsarextra/genuimodels

    """

    from pandas import DataFrame

    from genui.qsar.genuimodels.bases import DescriptorCalculator

    from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator
    from rdkit.Chem import Descriptors, MolFromSmiles

    class RDKitDescriptorsCalculator(DescriptorCalculator):
        group_name = 'RDKit_2D'

        def __call__(self, smiles) -> DataFrame:
            """
            Calculates 2D RDKit descriptors.

            Parameters
            ----------
            smiles : list
                A list of SMILES strings.

            Returns
            -------
            descriptors : DataFrame
                The matrix of calculated descriptors as `DataFrame`.
            """

            desc_list = [x[0] for x in Descriptors.descList]
            calc = MolecularDescriptorCalculator(desc_list)
            ret = []
            for smile in smiles:
                mol = MolFromSmiles(smile)
                descs = calc.CalcDescriptors(mol)
                ret.append(descs)

            return DataFrame(ret, columns=desc_list)

Note that you also have to give the new group of descriptors a name using the `DescriptorCalculator.group_name` class attribute. This is the name under
which this descriptor group appears in the REST API.

..  _dev-guide-qsar-metrics:

Adding Performance Metrics
--------------------------

GenUI already has a small collection of performance metrics for both
classification and regression tasks. However, it is very easy to
implement custom metrics. The process is similar to what we have
seen so far. You just need to create a new :code:`metrics.py`
module file in the :code:`qsarextra.genuimodels` package
and create subclasses of `ValidationMetric` inside.

We could again exploit scikit-learn to provide a simple
implementation of the F1 score:

..  code-block:: python

    """
    metrics.py in src/genui/qsar/extensions/qsarextra/genuimodels/
    """

    from sklearn import metrics

    from genui.models.genuimodels.bases import ValidationMetric, Algorithm


    class F1(ValidationMetric):
        """
        Implementation of the F1 score for classification accuracy.
        """

        name = "F1"
        description = "Compute the F1 score, also known as balanced F-score or F-measure."
        modes = [Algorithm.CLASSIFICATION]

        def __call__(self, true_vals, predicted_vals):
            """
            Implementation of the validation metric calculation.

            Note: Predicted values (`predicted_vals`) for classification models
            should be probabilities with which an item belongs
            to the class noted in `true_vals`. For regression, `predicted_vals`
            are simply the predicted values.

            Parameters
            ----------
            true_vals
                True prediction values from data.
            predicted_vals
                Predicted values from the model.
            Returns
            -------
            score :float
                A single number representing the model score according to this metric.

            """

            return metrics.f1_score(true_vals, self.probasToClasses(predicted_vals))

All you have to do is implement the :code:`__call__` method and give your new metric
a name, description and a list of modes you want this metric to be available for.