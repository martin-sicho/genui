import json
import os
import tarfile
import tempfile

from django.core.exceptions import ImproperlyConfigured
from django.urls import reverse
from qsprpred.models import SklearnModel
from rest_framework.test import APITestCase
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

from genui.compounds.extensions.chembl.tests import CompoundsMixIn
from genui.compounds.models import ActivityTypes, ActivityUnits
from genui.models.models import ModelPerformance, Algorithm, AlgorithmMode, ModelFile, \
    BasicValidationStrategy, MetricCurvePoint
from genui.qsar.models import QSARModel, ModelActivitySet
from .genuimodels import builders


class QSARModelInit(CompoundsMixIn):
    """
    Base class for initializing QSAR model tests.
    
    This class sets up the necessary data structures and methods
    for creating and testing QSAR models.
    """

    def setUp(self):
        super().setUp()
        self.project = self.createProject()
        self.molset = self.createMolSet(
            reverse('chemblSet-list'),
            {
                "targets": ["CHEMBL251"],
                "maxPerTarget": 30
            }
        )

    def createTestQSARModel(

            self,
            activitySet=None,
            activityType=None,
            mode=None,
            algorithm=None,
            parameters=None,
            embeddings=None,
            metrics=None,
            dataSplit=None,
            hyperParamOptStrategies=None,
            correct=True,
            validationStrategies=None,
    ):
        """
        Create a test QSAR model with specified parameters.

        Args:
            activitySet (ActivitySet, optional): The activity set to use.
            activityType (ActivityType, optional): The type of activity.
            mode (AlgorithmMode, optional): The mode of the algorithm.
            algorithm (Algorithm, optional): The algorithm to use.
            parameters (dict, optional): Algorithm parameters.
            embeddings (list, optional): List of embeddings.
            metrics (list, optional): List of performance metrics.

        Returns:
            QSARModel: The created QSAR model instance.
        """

        if not activitySet:
            activitySet = self.molset.activities.all()[0]
        if not activityType:
            activityType = ActivityTypes.objects.get(value="Ki_pChEMBL")
        if not mode:
            mode = AlgorithmMode.objects.get(name="classification")
        if not algorithm:
            algorithm = Algorithm.objects.get(name="QSPRPredScikitModel")
        if not parameters:
            parameters = {"alg": "RandomForestClassifier",
                          "parameters": json.dumps({"n_estimators": 150, })
                          }
        if not embeddings:
            embeddings = [{"name": "MorganFP", "arguments": {"radius": 2, "nBits": 2048}}, ]
        if not metrics:
            metrics = ["matthews_corrcoef", "roc_curve"]
        if not dataSplit:
            dataSplit = {"name": "qsprpred.data.sampling.splits.RandomSplit", "test_fraction": 0.2, "seed": 42}
        if not hyperParamOptStrategies:
            hyperParamOptStrategies = []

        post_data = {
            "name": "Test Model",
            "description": "test description",
            "project": self.project.id,
            "molset": self.molset.id,
            "trainingStrategy": {
                "algorithm": algorithm.id,
                "parameters": parameters,
                "mode": mode.id,
                "embeddings": embeddings,
                "activityThreshold": 6.5,
                "activitySet": activitySet.id,
                "activityType": activityType.id,
                "validationStrategies": [{
                    "resourcetype": "BasicValidationStrategy",
                    "dataSplit": dataSplit,
                    "cvFolds": 3,
                    "metrics": metrics
                }] if not validationStrategies else validationStrategies,
                "hyperParamOptStrategies": hyperParamOptStrategies
            }
        }
        create_url = reverse('model-list')
        response = self.client.post(create_url, data=post_data, format='json')
        print(json.dumps(response.data, indent=4))

        if correct:
            self.assertEqual(response.status_code, 201)
            return QSARModel.objects.get(pk=response.data["id"])
        else:
            return response

    def predictWithModel(self, model, to_predict):
        """
        Make predictions using the given model on the specified molecules.

        Args:
            model (QSARModel): The QSAR model to use for predictions.
            to_predict (MolSet): The set of molecules to predict.

        Returns:
            ModelActivitySet: The resulting set of predicted activities.
        """

        post_data = {
            "name": f"Predictions using {model.name}",
            "molecules": to_predict.id
        }
        create_url = reverse('model-predictions', args=[model.id])
        response = self.client.post(create_url, data=post_data, format='json')
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 201)

        instance = ModelActivitySet.objects.get(pk=response.data['id'])
        url = reverse('activitySet-activities', args=[instance.id])
        response = self.client.get(url)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.data['count'], to_predict.molecules.count())
        print(json.dumps(response.data, indent=4))

        return instance

    def uploadModel(self, filePath, algorithm, mode, predictionsType, predictionsUnits):
        """
        Upload a pre-trained model file and create a corresponding QSAR model.

        Args:
            filePath (str): Path to the model file.
            algorithm (Algorithm): The algorithm used in the model.
            mode (AlgorithmMode): The mode of the algorithm.
            predictionsType (str): Type of predictions the model makes.
            predictionsUnits (str): Units of the predictions.

        Returns:
            QSARModel: The created QSAR model instance.
        """

        create_url = reverse('model-list')
        post_data = {
            "name": "Test Model",
            "description": "test description",
            "project": self.project.id,
            "build": False,
            "predictionsType": predictionsType,
            "predictionsUnits": predictionsUnits,
            "trainingStrategy": {
                "algorithm": algorithm.id,
                "mode": mode.id,
            },
        }
        response = self.client.post(create_url, data=post_data, format='json')
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 201)
        instance = QSARModel.objects.get(pk=response.data["id"])
        self.assertFalse(instance.modelFile)

        url = reverse('qsar-model-files-list', args=[instance.id])
        response = self.client.post(
            url,
            data={
                "file": open(filePath, "rb"),
                "kind": ModelFile.MAIN,
            },
            format='multipart'
        )
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 201)

        url = reverse('model-detail', args=[instance.id])
        response_other = self.client.get(url)
        self.assertEqual(response.data['file'].split('/')[-1], response_other.data['modelFile']['file'].split('/')[-1])

        return instance


class ModelInitTestCase(QSARModelInit, APITestCase):
    """Test case for QSAR model initialization and basic functionality."""

    def test_create_view_classification(self):
        """Test the creation and basic functionality of a classification QSAR model."""
        model = self.createTestQSARModel()

        path = model.modelFile.path
        temp_dir = tempfile.TemporaryDirectory()
        with tarfile.open(path, 'r:*') as tar_ref:
            tar_ref.extractall(temp_dir.name)

        model_alg_major = "QSPRPredScikitModel"
        model_alg_minor = "RandomForestClassifier"
        model_name = f"{model_alg_major}_{model_alg_minor}"
        model = SklearnModel.fromFile(os.path.join(temp_dir.name, model_name, model_name + "_meta.json"))
        self.assertTrue(isinstance(model.estimator, RandomForestClassifier))

        # get the model via api
        response = self.client.get(reverse('model-list'))
        self.assertEqual(response.status_code, 200)
        print(json.dumps(response.data[0], indent=4))

        # create predictions with the model
        model = QSARModel.objects.get(pk=response.data[0]['id'])
        self.predictWithModel(model, self.molset)

        # make sure the delete cascades fine and the file gets deleted too
        self.project.delete()
        self.assertTrue(ModelPerformance.objects.count() == 0)
        self.assertTrue(not os.path.exists(path))

    def test_create_view_from_file_classification(self):
        """Test creating a classification QSAR model from a pre-trained model file."""
        instance_first = self.createTestQSARModel()
        self.assertEqual(instance_first.predictionsType, ActivityTypes.objects.get(value="Active Probability"))
        self.assertEqual(instance_first.predictionsUnits, None)
        instance = self.uploadModel(
            instance_first.modelFile.path,
            instance_first.trainingStrategy.algorithm,
            instance_first.trainingStrategy.mode,
            instance_first.predictionsType.value,
            instance_first.predictionsUnits.value if instance_first.predictionsUnits else None,
        )

        builder = builders.BasicQSARModelBuilder(instance)
        self.assertRaisesMessage(ImproperlyConfigured, "You cannot build a QSAR model without validation strategies.",
                                 builder.build)
        print(builder.predictMols(["CC", "CCO"]))

        activity_set = self.predictWithModel(instance, self.molset)
        for activity in activity_set.activities.all():
            self.assertEqual(activity.type, instance_first.predictionsType)
            self.assertEqual(activity.units, instance_first.predictionsUnits)

    def test_create_view_regression(self):
        """Test the creation and basic functionality of a regression QSAR model."""
        model = self.createTestQSARModel(
            mode=AlgorithmMode.objects.get(name="regression"),
            parameters={"alg": "RandomForestRegressor",
                        "parameters": json.dumps({"n_estimators": 150, })
                        },
            metrics=["r2", "neg_mean_squared_error"],
            activityType=ActivityTypes.objects.get(value="Ki")
        )
        self.assertEqual(model.predictionsType, ActivityTypes.objects.get(value="Ki"))
        self.assertEqual(model.predictionsUnits, ActivityUnits.objects.get(value="nM"))
        path = model.modelFile.path
        temp_dir = tempfile.TemporaryDirectory()
        with tarfile.open(path, 'r:*') as tar_ref:
            tar_ref.extractall(temp_dir.name)

        model_alg_major = "QSPRPredScikitModel"
        model_alg_minor = "RandomForestRegressor"
        model_name = f"{model_alg_major}_{model_alg_minor}"
        model_ = SklearnModel.fromFile(os.path.join(temp_dir.name, model_name, model_name + "_meta.json"))
        self.assertTrue(isinstance(model_.estimator, RandomForestRegressor))
        activity_set_orig = self.predictWithModel(model, self.molset)

        # try to upload it as a file and use that model for predictions
        model_from_file = self.uploadModel(
            model.modelFile.path,
            model.trainingStrategy.algorithm,
            model.trainingStrategy.mode,
            model.predictionsType.value,
            model.predictionsUnits.value if model.predictionsUnits else None,
        )
        builder = builders.BasicQSARModelBuilder(model_from_file)
        print(builder.predictMols(["CC", "CCO"]))

        activity_set = self.predictWithModel(model_from_file, self.molset)
        for activity_uploaded, activity_orig in zip(activity_set.activities.all(), activity_set_orig.activities.all()):
            self.assertEqual(activity_uploaded.type, model.predictionsType)
            self.assertEqual(activity_uploaded.units, model.predictionsUnits)
            self.assertEqual(activity_uploaded.type, activity_orig.type)
            self.assertEqual(activity_uploaded.units, activity_orig.units)
            self.assertEqual(activity_uploaded.value, activity_orig.value)

    def test_training_strategy_has_validation_strategies(self):
        """Test that the training strategy of a QSAR model has validation strategies."""
        # Create a QSAR model
        model = self.createTestQSARModel()

        # Check if the training strategy has validation strategies
        self.assertTrue(model.trainingStrategy.validationStrategies.exists())
        self.assertEqual(model.trainingStrategy.validationStrategies.count(), 1)

    def test_multiple_validation_strategies(self):
        """Test adding multiple validation strategies to a QSAR model."""
        # Create initial QSAR model with one validation strategy
        model = self.createTestQSARModel()

        # Add a second validation strategy
        ds = {"name": "qsprpred.data.sampling.splits.RandomSplit",
              "test_fraction": 0.2,
              "seed": 42,
              }
        second_strategy = BasicValidationStrategy.objects.create(
            trainingStrategy=model.trainingStrategy,
            cvFolds=5,
            dataSplit=ds,
        )
        second_strategy.metrics = ["r2", "mean_squared_error"]
        model.trainingStrategy.validationStrategies.add(second_strategy)

        # Check if the training strategy has multiple validation strategies
        self.assertEqual(model.trainingStrategy.validationStrategies.count(), 2)

        # Verify that the validation strategies are different
        validation_strategies = list(model.trainingStrategy.validationStrategies.all())
        self.assertNotEqual(validation_strategies[0].cvFolds, validation_strategies[1].cvFolds)
        self.assertNotEqual(set(validation_strategies[0].metrics), set(validation_strategies[1].metrics))

    def test_default_validation_strategy_parameters(self):
        """
        Test that the default validation strategy parameters are set correctly.
                
        The default validation strategy should have the following parameters:
        - cvFolds: 3
        - validSetSize: 0.2
        - metrics: matthews_corrcoef
        """
        model = self.createTestQSARModel()
        validation_strategy = model.trainingStrategy.validationStrategies.first()
        self.assertEqual(validation_strategy.cvFolds, 3)
        self.assertEqual(validation_strategy.dataSplit["test_fraction"], 0.2)
        self.assertEqual(set(validation_strategy.metrics), set(["matthews_corrcoef", "roc_curve"]))

    def test_update_validation_strategy(self):
        """Test that the validation strategy parameters can be updated"""
        model = self.createTestQSARModel()
        validation_strategy = model.trainingStrategy.validationStrategies.first()
        validation_strategy.cvFolds = 10
        validation_strategy.save()
        updated_strategy = BasicValidationStrategy.objects.get(id=validation_strategy.id)
        self.assertEqual(updated_strategy.cvFolds, 10)

    def test_remove_validation_strategy(self):
        """Test removing a validation strategy from a QSAR model."""
        model = self.createTestQSARModel()
        validation_strategy = model.trainingStrategy.validationStrategies.first()
        validation_strategy.delete()
        self.assertFalse(model.trainingStrategy.validationStrategies.exists())

    def test_different_models_different_validation_strategies(self):
        """Test that different QSAR models have different validation strategies."""
        model1 = self.createTestQSARModel()
        model2 = self.createTestQSARModel()
        strategy1 = model1.trainingStrategy.validationStrategies.first()
        strategy2 = model2.trainingStrategy.validationStrategies.first()
        self.assertNotEqual(strategy1, strategy2)

    def test_performance_metrics_associated_with_validation_strategies(self):
        """Test associating performance metrics with validation strategies."""
        model = self.createTestQSARModel()
        validation_strategy = model.trainingStrategy.validationStrategies.first()
        metrics = ["r2", "mean_squared_error"]
        validation_strategy.metrics = metrics
        validation_strategy.save()
        self.assertEqual(set(validation_strategy.metrics), set(metrics))

    def test_change_data_split(self):
        model = self.createTestQSARModel()
        validation_strategy = model.trainingStrategy.validationStrategies.first()
        new_split = {"name": "qsprpred.data.sampling.splits.RandomSplit",
                     "test_fraction": 0.3,
                     "seed": 69,
                     }
        validation_strategy.dataSplit = new_split
        validation_strategy.save()
        self.assertEqual(model.trainingStrategy.validationStrategies.first().dataSplit["test_fraction"], 0.3)
        self.assertEqual(model.trainingStrategy.validationStrategies.first().dataSplit["seed"], 69)

    def test_qsprpred_fingerprints(self):
        fingerprints = [
            {"name": "MorganFP", "arguments": {"radius": 2, "nBits": 2048, }},
            {"name": "RDKitMACCSFP", "arguments": {}},
            {"name": "MaccsFP", "arguments": {"nBits": 167}},
            {"name": "AvalonFP", "arguments": {"nBits": 1024}},
            {"name": "TopologicalFP", "arguments": {"nBits": 2048}},
            {"name": "AtomPairFP", "arguments": {"nBits": 2048}},
            {"name": "RDKitFP", "arguments": {"minPath": 1, "maxPath": 7, "nBits": 2048}},
            {"name": "PatternFP", "arguments": {"nBits": 2048}},
            {"name": "LayeredFP", "arguments": {"minPath": 1, "maxPath": 7, "nBits": 2048}},
        ]
        model = self.createTestQSARModel(embeddings=fingerprints)

        response = self.client.get(reverse('model-list'))
        self.assertEqual(response.status_code, 200)
        print(json.dumps(response.data[0], indent=4))

        # create predictions with the model
        model = QSARModel.objects.get(pk=response.data[0]['id'])
        self.predictWithModel(model, self.molset)

    def test_qsprpred_embedding_set(self):
        sets = [{"name": "DrugExPhyschem", "arguments": {"physchem_props": ["Aliphatic"]}},
                {"name": "RDKitDescs", "arguments": {"rdkit_descriptors": ["NumSaturatedCarbocycles"]}}, ]
        model = self.createTestQSARModel(embeddings=sets)

        response = self.client.get(reverse('model-list'))
        self.assertEqual(response.status_code, 200)
        print(json.dumps(response.data[0], indent=4))

        # create predictions with the model
        model = QSARModel.objects.get(pk=response.data[0]['id'])
        self.predictWithModel(model, self.molset)

    def test_change_algorithm(self):
        parameters = {"alg": "KNeighborsClassifier", "parameters": json.dumps({"n_neighbors": 5})}
        model = self.createTestQSARModel(parameters=parameters)

    def test_incorrect_algorithm_parameters_and_combinations(self):
        response = self.createTestQSARModel(parameters={"alg": "KNeighborsClassifier",
                                                        "parameters": json.dumps({"n_estimators": 5})},
                                            correct=False)
        self.assertEqual(response.status_code, 400)
        self.assertEqual('Parameter n_estimators is not valid for the selected algorithm KNeighborsClassifier.',
                         response.json()[0])

        response = self.createTestQSARModel(mode=AlgorithmMode.objects.get(name="regression"), correct=False)
        self.assertEqual(response.status_code, 400)
        self.assertEqual('You cannot use a classifier algorithm for a regression model.',
                         response.json()[0])

    def test_upload_external_model_file(self):
        from qsprpred.data import QSPRDataset
        from qsprpred.data import RandomSplit
        from qsprpred.data.descriptors.fingerprints import MorganFP
        import pandas as pd

        self.molset = self.createMolSet(
            reverse('chemblSet-list'),
            {
                "targets": ["CHEMBL251"],
                "maxPerTarget": 50
            }
        )
        temp_dir_model = tempfile.TemporaryDirectory()

        activity_set = self.molset.activities.get()
        compounds, activities, units = activity_set.cleanForModelling(ActivityTypes.objects.get(value="Ki_pChEMBL").id)
        compounds = [x.canonicalSMILES for x in compounds]

        dataset = QSPRDataset("TestDataset",
                              [{"name": "activity", "task": "SINGLECLASS", "th": [6.5]}],
                              pd.DataFrame({'SMILES': compounds, 'activity': activities}))

        rand_split = RandomSplit(test_fraction=0.2, dataset=dataset)
        dataset.prepareDataset(
            split=rand_split,
            feature_calculators=[MorganFP(radius=3, nBits=2048)],
        )
        alg = SklearnModel(temp_dir_model.name, RandomForestClassifier,
                           "TestModel", parameters={"n_estimators": 150})
        alg.fitDataset(dataset)
        alg.save(True)
        model_path = os.path.join(temp_dir_model.name, alg.name + ".tar.gz")
        with tarfile.open(model_path, "w:gz") as tar:
            tar.add(os.path.join(temp_dir_model.name, alg.name), arcname=alg.name)
        instance = self.uploadModel(
            model_path,
            Algorithm.objects.get_or_create(name="QSPRPredScikitModel")[0],
            AlgorithmMode.objects.get(name="classification"),
            'Active Probability',
            None,
        )
        self.predictWithModel(instance, self.molset)

    def test_hyperparameter_optimization(self):
        self.createTestQSARModel(hyperParamOptStrategies=[{
            "resourcetype": "GridSearchOptimization",
            "searchSpace": [{"name": "n_estimators", "type": "range", "value": [100, 150, 10]},
                            {"name": "max_depth", "type": "range", "value": [5, 10]},
                            {"name": "criterion", "type": "sequence", "value": ["gini", "entropy", "log_loss"]}],
            "metric": "f1",
            "scoreAggregation": "mean"
        }], )
        self.createTestQSARModel(hyperParamOptStrategies=[{
            "resourcetype": "OptunaOptimization",
            "searchSpace": [{"name": "n_estimators", "type": "int", "value": [100, 250]},
                            {"name": "criterion", "type": "categorical", "value": ["gini", "entropy", "log_loss"]}],
            "metric": "f1",
            "scoreAggregation": "mean",
            "nTrials": 10,
        }])

    def test_multiple_validation_strategies_truly(self):
        validation_strategies = [
            {
                "resourcetype": "BasicValidationStrategy",
                "dataSplit": {"name": "qsprpred.data.sampling.splits.RandomSplit",
                              "test_fraction": 0.2,
                              "seed": 42,
                              },
                "cvFolds": 3,
                "metrics": ["matthews_corrcoef"]
            },
            {
                "resourcetype": "BasicValidationStrategy",
                "dataSplit": {"name": "qsprpred.data.sampling.splits.RandomSplit",
                              "test_fraction": 0.25,
                              "seed": 111,
                              },
                "cvFolds": 5,
                "metrics": ["matthews_corrcoef", "accuracy"]
            }
        ]
        model = self.createTestQSARModel(validationStrategies=validation_strategies)

    def test_all_parameters_default_value(self):  # GUI creates a problem with default values, Here's an example
        model = self.createTestQSARModel(
            parameters={"alg": "RandomForestClassifier",
                        "parameters": json.dumps({"n_estimators": 100,
                                                  "criterion": "gini",
                                                  "max_depth": 1,
                                                  "min_samples_split": 2,
                                                  "min_weight_fraction_leaf": 0.0,
                                                  "max_features": "sqrt",
                                                  "min_impurity_decrease": 0,
                                                  "bootstrap": True,
                                                  "oob_score": True,
                                                  "warm_start": True,
                                                  "class_weight": "balanced",
                                                  "ccp_alpha": 0.0,
                                                  "max_samples": 0.001})
                        }, )

    def test_all_metric_curves(self):
        model = self.createTestQSARModel(metrics=["roc_curve", "precision_recall_curve", "det_curve"])
        print(MetricCurvePoint.objects.filter(model=model))

    def test_compound_splits(self):
        metrics = ["roc_curve", "accuracy"]
        bootstrap_split = {
            "name": "qsprpred.data.sampling.splits.BootstrapSplit",
            "split": {"name": "qsprpred.data.sampling.splits.RandomSplit",
                      "test_fraction": 0.2,
                      "seed": 42,
                      },
            "n_bootstraps": 3
        }
        gbmt_data_split = {"name": "qsprpred.data.sampling.splits.GBMTDataSplit",
                           "clustering":
                               {"name": "qsprpred.data.chem.clustering.RandomClusters",
                                }
                           }
        scaffold_split = {"name": "qsprpred.data.sampling.splits.ScaffoldSplit",
                          "scaffold":
                              {"name": "qsprpred.data.chem.scaffolds.BemisMurcko"
                               },
                          "n_folds": 2
                          }
        cluster_split = {"name": "qsprpred.data.sampling.splits.ClusterSplit",
                         "clustering": {
                             "name": "qsprpred.data.chem.clustering.FPSimilarityMaxMinClusters",
                             "fp_calculator": {
                                 "name": "qsprpred.data.descriptors.fingerprints.AvalonFP",
                                 "nBits": 512
                             }
                         },
                         "seed": 123
                         }
        validation_strategies = [
            {
                "resourcetype": "BasicValidationStrategy",
                "dataSplit": bootstrap_split,
                "cvFolds": 3,
                "metrics": metrics
            },
            {
                "resourcetype": "BasicValidationStrategy",
                "dataSplit": gbmt_data_split,
                "cvFolds": 3,
                "metrics": metrics
            },
            {
                "resourcetype": "BasicValidationStrategy",
                "dataSplit": scaffold_split,
                "cvFolds": 3,
                "metrics": metrics
            },
            {
                "resourcetype": "BasicValidationStrategy",
                "dataSplit": cluster_split,
                "cvFolds": 3,
                "metrics": metrics
            }
        ]
        model = self.createTestQSARModel(validationStrategies=validation_strategies, )

    def test_sklearn_split(self):
        data_split = {"name": "sklearn.model_selection.ShuffleSplit",
                      "test_size": 0.2,
                      }
        model = self.createTestQSARModel(dataSplit=data_split)

    def test_xgboost(self):
        model = self.createTestQSARModel(
            parameters={"alg": "xgboost.XGBClassifier",
                        "parameters": json.dumps({})
                        },)

        model = self.createTestQSARModel(
            mode=AlgorithmMode.objects.get(name="regression"),
            parameters={"alg": "xgboost.XGBRegressor",
                        "parameters": json.dumps({})
                        },
            metrics=["r2", "neg_mean_squared_error"],
            activityType=ActivityTypes.objects.get(value="Ki")
        )
