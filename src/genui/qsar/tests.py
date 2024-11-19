import json
import os

import joblib
from django.core.exceptions import ImproperlyConfigured
from rest_framework.test import APITestCase
from django.urls import reverse
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

from genui.compounds.models import ActivityTypes, ActivityUnits
from genui.compounds.extensions.chembl.tests import CompoundsMixIn
from genui.qsar.models import QSARModel, DescriptorGroup, ModelActivitySet, TrainingStrategy
from genui.models.models import ModelPerformance, Algorithm, AlgorithmMode, ModelFile, ModelPerformanceMetric, BasicValidationStrategy
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
                "maxPerTarget" : 30
            }
        )

    def createTestQSARModel(
        
            self,
            activitySet=None,
            activityType=None,
            mode=None,
            algorithm=None,
            parameters=None,
            descriptors=None,
            metrics=None
    ):
        """
        Create a test QSAR model with specified parameters.

        Args:
            activitySet (ActivitySet, optional): The activity set to use.
            activityType (ActivityType, optional): The type of activity.
            mode (AlgorithmMode, optional): The mode of the algorithm.
            algorithm (Algorithm, optional): The algorithm to use.
            parameters (dict, optional): Algorithm parameters.
            descriptors (list, optional): List of descriptors.
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
            algorithm = Algorithm.objects.get(name="RandomForest")
        if not parameters:
            parameters = {
                "n_estimators": 150
            }
        if not descriptors:
            descriptors = [DescriptorGroup.objects.get(name="MORGANFP")]
        if not metrics:
            metrics = [
                ModelPerformanceMetric.objects.get(name="MCC"),
                ModelPerformanceMetric.objects.get(name="ROC"),
            ]

        post_data = {
            "name": "Test Model",
            "description": "test description",
            "project": self.project.id,
            "molset": self.molset.id,
            "trainingStrategy": {
                "algorithm": algorithm.id,
                "parameters": parameters,
                "mode": mode.id,
                "descriptors": [
                    x.id for x in descriptors
                ],
                "activityThreshold": 6.5,
                "activitySet": activitySet.id,
                "activityType": activityType.id,
                "validationStrategies": [{
                    "resourcetype": "BasicValidationStrategy",
                    "cvFolds": 3,
                    "validSetSize": 0.2,
                    "metrics": [
                        x.id for x in metrics
                    ]
                }]
            }
        }
        create_url = reverse('model-list')
        response = self.client.post(create_url, data=post_data, format='json')
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 201)

        return QSARModel.objects.get(pk=response.data["id"])

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

    def uploadModel(self, filePath, algorithm, mode, descriptors, predictionsType, predictionsUnits):
        """
        Upload a pre-trained model file and create a corresponding QSAR model.

        Args:
            filePath (str): Path to the model file.
            algorithm (Algorithm): The algorithm used in the model.
            mode (AlgorithmMode): The mode of the algorithm.
            descriptors (list): List of descriptors used in the model.
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
            "build" : False,
            "predictionsType": predictionsType,
            "predictionsUnits": predictionsUnits,
            "trainingStrategy": {
                "algorithm": algorithm.id,
                "mode": mode.id,
                "descriptors": [
                  x.id for x in descriptors
                ]
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
                "file" : open(filePath, "rb"),
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
        model = joblib.load(model.modelFile.path)
        self.assertTrue(isinstance(model, RandomForestClassifier))

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
            [DescriptorGroup.objects.get(name='MORGANFP')],
            instance_first.predictionsType.value,
            instance_first.predictionsUnits.value if instance_first.predictionsUnits else None
        )

        builder = builders.BasicQSARModelBuilder(instance)
        self.assertRaisesMessage(ImproperlyConfigured, "You cannot build a QSAR model without validation strategies.", builder.build)
        builder.calculateDescriptors(["CC", "CCO"])
        print(builder.predict())

        activity_set = self.predictWithModel(instance, self.molset)
        for activity in activity_set.activities.all():
            self.assertEqual(activity.type, instance_first.predictionsType)
            self.assertEqual(activity.units, instance_first.predictionsUnits)

    def test_create_view_regression(self):
        """Test the creation and basic functionality of a regression QSAR model."""
        model = self.createTestQSARModel(
            mode=AlgorithmMode.objects.get(name="regression"),
            metrics=ModelPerformanceMetric.objects.filter(name__in=("R2", "MSE")),
            activityType=ActivityTypes.objects.get(value="Ki")
        )
        self.assertEqual(model.predictionsType, ActivityTypes.objects.get(value="Ki"))
        self.assertEqual(model.predictionsUnits, ActivityUnits.objects.get(value="nM"))
        self.assertTrue(isinstance(joblib.load(model.modelFile.path), RandomForestRegressor))
        activity_set_orig = self.predictWithModel(model, self.molset)

        # try to upload it as a file and use that model for predictions
        model_from_file = self.uploadModel(
            model.modelFile.path,
            model.trainingStrategy.algorithm,
            model.trainingStrategy.mode,
            [DescriptorGroup.objects.get(name='MORGANFP')],
            model.predictionsType.value,
            model.predictionsUnits.value if model.predictionsUnits else None
        )
        builder = builders.BasicQSARModelBuilder(model_from_file)
        builder.calculateDescriptors(["CC", "CCO"])
        print(builder.predict())
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
        second_strategy = BasicValidationStrategy.objects.create(
            trainingStrategy=model.trainingStrategy,
            cvFolds=5,
            validSetSize=0.2
        )
        second_strategy.metrics.set(ModelPerformanceMetric.objects.filter(name__in=["R2", "MSE"]))
        model.trainingStrategy.validationStrategies.add(second_strategy)
        
        # Check if the training strategy has multiple validation strategies
        self.assertEqual(model.trainingStrategy.validationStrategies.count(), 2)
        
        # Verify that the validation strategies are different
        validation_strategies = list(model.trainingStrategy.validationStrategies.all())
        self.assertNotEqual(validation_strategies[0].cvFolds, validation_strategies[1].cvFolds)
        self.assertNotEqual(set(validation_strategies[0].metrics.all()), set(validation_strategies[1].metrics.all()))

    def test_default_validation_strategy_parameters(self):
        """
        Test that the default validation strategy parameters are set correctly.
                
        The default validation strategy should have the following parameters:
        - cvFolds: 3
        - validSetSize: 0.2
        - metrics: MCC, ROC
        """        
        model = self.createTestQSARModel()
        validation_strategy = model.trainingStrategy.validationStrategies.first()
        self.assertEqual(validation_strategy.cvFolds, 3)
        self.assertEqual(validation_strategy.validSetSize, 0.2)
        self.assertEqual(set(validation_strategy.metrics.all()), set(ModelPerformanceMetric.objects.filter(name__in=["MCC", "ROC"]))
        )

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
        metrics = ModelPerformanceMetric.objects.filter(name__in=["R2", "MSE"])
        validation_strategy.metrics.set(metrics)
        validation_strategy.save()
        self.assertEqual(set(validation_strategy.metrics.all()), set(metrics))
        