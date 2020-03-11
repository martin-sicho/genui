import json
import os

import joblib
from django.core.exceptions import ImproperlyConfigured
from rest_framework.test import APITestCase, APITransactionTestCase
from django.urls import reverse
from sklearn.ensemble import RandomForestClassifier

from compounds.initializers.chembl import ChEMBLSetInitializer
from compounds.models import ChEMBLCompounds
from modelling.apps import ModellingConfig
from projects.models import Project
from qsar.models import QSARModel, DescriptorGroup, ModelActivitySet
from modelling.models import ModelPerformance, Algorithm, AlgorithmMode, ModelFile, ModelPerformanceMetric
from .core import builders


class InitMixIn:

    def setUp(self):
        from qsar.apps import QsarConfig
        ModellingConfig.ready('dummy', True)
        QsarConfig.ready('dummy', True)
        self.project = Project.objects.create(**{
            "name" : "Test Project"
            , "description" : "Test Description"
        })
        self.molset = ChEMBLCompounds.objects.create(**{
            "name": "Test ChEMBL Data Set",
            "description": "Some description...",
            "project": self.project
        })
        initializer = ChEMBLSetInitializer(
            self.molset
            , targets=["CHEMBL251"]
            , max_per_target=50
        )
        initializer.populateInstance()

    def tearDown(self) -> None:
        if self.project.id:
            self.project.delete()

    def createTestModel(self):
        post_data = {
          "name": "Test Model",
          "description": "test description",
          "project": self.project.id,
          "trainingStrategy": {
            "algorithm": Algorithm.objects.get(name="RandomForest").id,
            "parameters": {
              "n_estimators": 150
            },
            "mode": AlgorithmMode.objects.get(name="classification").id,
            "descriptors": [
              DescriptorGroup.objects.get(name="MORGANFP").id
            ],
            "activityThreshold": 6.5
          },
          "validationStrategy": {
            "cvFolds": 3,
            "validSetSize": 0.2,
            "metrics": [
              ModelPerformanceMetric.objects.get(name="MCC").id
            ]
          },
          "molset": self.molset.id
        }
        create_url = reverse('model-list')
        response = self.client.post(create_url, data=post_data, format='json')
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 201)

        instance = QSARModel.objects.get(pk=response.data["id"])
        builder_class = 'BasicQSARModelBuilder'
        builder_class = getattr(builders, builder_class)
        builder = builder_class(instance)
        builder.build()

        return instance

class ModelInitTestCase(InitMixIn, APITestCase):

    def test_create_view(self):
        instance = self.createTestModel()

        path = instance.modelFile.path
        model = joblib.load(instance.modelFile.path)
        self.assertTrue(isinstance(model, RandomForestClassifier))

        # get the model via api
        response = self.client.get(reverse('model-list'))
        self.assertEqual(response.status_code, 200)
        print(json.dumps(response.data[0], indent=4))

        # create predictions with the model
        model = QSARModel.objects.get(pk=response.data[0]['id'])
        post_data = {
            "name": f"Predictions using {model.name}",
            "molecules": self.molset.id
        }
        create_url = reverse('model-predictions', args=[model.id])
        response = self.client.post(create_url, data=post_data, format='json')
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 201)

        instance = ModelActivitySet.objects.get(pk=response.data['id'])
        model = QSARModel.objects.get(pk=instance.model.id)
        builder_class = getattr(builders, model.builder.name)
        builder = builder_class(
            model
        )
        builder.populateActivitySet(instance)

        url = reverse('activitySet-activities', args=[response.data['id']])
        response = self.client.get(url)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.data['count'], self.molset.molecules.count())
        print(json.dumps(response.data, indent=4))

        # make sure the delete cascades fine and the file gets deleted too
        self.project.delete()
        self.assertTrue(ModelPerformance.objects.count() == 0)
        self.assertTrue(not os.path.exists(path))

    def test_create_view_from_file(self):
        instance_first = self.createTestModel()
        upload_file = open(instance_first.modelFile.path, "rb")

        create_url = reverse('model-list')
        post_data = {
            "name": "Test Model",
            "description": "test description",
            "project": self.project.id,
            "build" : False,
            "trainingStrategy": {
                "algorithm": Algorithm.objects.get(name="RandomForest").id,
                "mode": AlgorithmMode.objects.get(name="classification").id,
                "descriptors": [
                  DescriptorGroup.objects.get(name="MORGANFP").id
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
                "file" : upload_file,
                "kind": ModelFile.MAIN,
            },
            format='multipart'
        )
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 201)

        response_other = self.client.get(reverse('model-list'), args=instance.id)
        self.assertEqual(response.data['file'].split('/')[-1], response_other.data[1]['modelFile']['file'].split('/')[-1])

        builder_class = 'BasicQSARModelBuilder'
        builder_class = getattr(builders, builder_class)
        builder = builder_class(instance)
        self.assertRaisesMessage(ImproperlyConfigured, "You cannot build a QSAR model with a missing validation strategy.", builder.build)
        builder.calculateDescriptors(["CC", "CCO"])
        print(builder.predict())