import json
import os

import joblib
from django.core.exceptions import ImproperlyConfigured
from rest_framework.test import APITestCase
from django.urls import reverse
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

from genui.compounds.models import ActivityTypes, ActivityUnits
from genui.compounds.extensions.chembl.tests import CompoundsMixIn
from genui.qsar.models import QSARModel, DescriptorGroup, ModelActivitySet
from genui.modelling.models import ModelPerformance, Algorithm, AlgorithmMode, ModelFile, ModelPerformanceMetric
from .genuimodels import builders


class QSARModelInit(CompoundsMixIn):

    def setUp(self):
        super().setUp()
        self.project = self.createProject()
        self.molset = self.createMolSet(
            reverse('chemblSet-list'),
            {
                "targets": ["CHEMBL251"],
                "maxPerTarget" : 100
            }
        )

    def tearDown(self) -> None:
        if self.project.id:
            self.project.delete()

    def createTestQSARModel(
            self,
            activitySet=None,
            activityType=None,
            mode=None,
            algorithm=None,
            descriptors=None,
            metrics=None
    ):
        if not activitySet:
            activitySet = self.molset.activities.all()[0]
        if not activityType:
            activityType = ActivityTypes.objects.get(value="Ki_pChEMBL")
        if not mode:
            mode = AlgorithmMode.objects.get(name="classification")
        if not algorithm:
            algorithm = Algorithm.objects.get(name="RandomForest")
        if not descriptors:
            descriptors = [DescriptorGroup.objects.get(name="MORGANFP")]
        if not metrics:
            metrics = [
                ModelPerformanceMetric.objects.get(name="MCC")
            ]

        post_data = {
            "name": "Test Model",
            "description": "test description",
            "project": self.project.id,
            "molset": self.molset.id,
            "trainingStrategy": {
                "algorithm": algorithm.id,
                "parameters": {
                    "n_estimators": 150
                },
                "mode": mode.id,
                "descriptors": [
                    x.id for x in descriptors
                ],
                "activityThreshold": 6.5,
                "activitySet": activitySet.id,
                "activityType": activityType.id
            },
            "validationStrategy": {
                "cvFolds": 3,
                "validSetSize": 0.2,
                "metrics": [
                    x.id for x in metrics
                ]
            }
        }
        create_url = reverse('model-list')
        response = self.client.post(create_url, data=post_data, format='json')
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 201)

        return QSARModel.objects.get(pk=response.data["id"])

    def predictWithModel(self, model, to_predict):
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

        response_other = self.client.get(reverse('model-list'), args=instance.id)
        self.assertEqual(response.data['file'].split('/')[-1], response_other.data[1]['modelFile']['file'].split('/')[-1])

        return instance

class ModelInitTestCase(QSARModelInit, APITestCase):

    def test_create_view_classification(self):
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
        instance_first = self.createTestQSARModel()
        self.assertEquals(instance_first.predictionsType, ActivityTypes.objects.get(value="Active Probability"))
        self.assertEquals(instance_first.predictionsUnits, None)
        instance = self.uploadModel(
            instance_first.modelFile.path,
            instance_first.trainingStrategy.algorithm,
            instance_first.trainingStrategy.mode,
            [DescriptorGroup.objects.get(name='MORGANFP')],
            instance_first.predictionsType.value,
            instance_first.predictionsUnits.value if instance_first.predictionsUnits else None
        )

        builder = builders.BasicQSARModelBuilder(instance)
        self.assertRaisesMessage(ImproperlyConfigured, "You cannot build a QSAR model with a missing validation strategy.", builder.build)
        builder.calculateDescriptors(["CC", "CCO"])
        print(builder.predict())

        activity_set = self.predictWithModel(instance, self.molset)
        for activity in activity_set.activities.all():
            self.assertEquals(activity.type, instance_first.predictionsType)
            self.assertEquals(activity.units, instance_first.predictionsUnits)

    def test_create_view_regression(self):
        model = self.createTestQSARModel(
            mode=AlgorithmMode.objects.get(name="regression"),
            metrics=ModelPerformanceMetric.objects.filter(name__in=("R2", "MSE")),
            activityType=ActivityTypes.objects.get(value="Ki")
        )
        self.assertEquals(model.predictionsType, ActivityTypes.objects.get(value="Ki"))
        self.assertEquals(model.predictionsUnits, ActivityUnits.objects.get(value="nM"))
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
            self.assertEquals(activity_uploaded.type, model.predictionsType)
            self.assertEquals(activity_uploaded.units, model.predictionsUnits)
            self.assertEquals(activity_uploaded.type, activity_orig.type)
            self.assertEquals(activity_uploaded.units, activity_orig.units)
            self.assertEquals(activity_uploaded.value, activity_orig.value)