import json
import os

import joblib
from rest_framework.test import APITestCase
from django.urls import reverse
from sklearn.ensemble import RandomForestClassifier

from compounds.initializers.chembl import ChEMBLSetInitializer
from compounds.models import ChEMBLCompounds
from modelling.apps import ModellingConfig
from projects.models import Project
from qsar.models import QSARModel
from modelling.models import ModelPerformance, Algorithm, AlgorithmMode
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
        self.post_data = {
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
              1
            ],
            "activityThreshold": 6.5
          },
          "validationStrategy": {
            "cvFolds": 10,
            "validSetSize": 0.2,
            "metrics": [
              1
            ]
          },
          "molset": self.molset.id
        }

    def tearDown(self) -> None:
        if self.project.id:
            self.project.delete()

    def createTestModel(self):
        create_url = reverse('model-list')
        response = self.client.post(create_url, data=self.post_data, format='json')
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
        model = joblib.load(instance.modelFile)
        self.assertTrue(isinstance(model, RandomForestClassifier))

        # get the model via api
        response = self.client.get(reverse('model-list'))
        print(json.dumps(response.data[0], indent=4))

        # make sure the delete cascades fine and the file gets deleted too
        self.project.delete()
        self.assertTrue(ModelPerformance.objects.count() == 0)
        self.assertTrue(not os.path.exists(path))


