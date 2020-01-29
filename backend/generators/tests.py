import json

from django.urls import reverse
from rest_framework.test import APITestCase

from generators.core import builders
from generators.models import DrugExNet
from modelling.models import Algorithm, AlgorithmMode
from qsar.tests import InitMixIn

class DrugExGeneratorInitTestCase(InitMixIn, APITestCase):

    def setUp(self):
        super().setUp()
        from generators.apps import GeneratorsConfig
        GeneratorsConfig.ready('dummy')
        # self.environ = self.createTestModel()

    def test_create_drugexnet_view(self):
        create_url = reverse("drugex_net-list")
        post_data = {
          "name": "Test DrugEx Network",
          "description": "test description",
          "project": self.project.id,
          "trainingStrategy": {
            "algorithm": Algorithm.objects.get(name="DrugExNetwork").id,
            "mode": AlgorithmMode.objects.get(name="generator").id,
            "parameters": {
                "nEpochs": 5,
                "monitorFrequency" : 10
            },
          },
          "validationStrategy": {
            "cvFolds": 10,
            "validSetSize": 0,
            "metrics": [
              1
            ]
          },
          "molset": self.molset.id
        }
        response = self.client.post(create_url, data=post_data, format='json')
        self.assertEqual(response.status_code, 201)
        print(json.dumps(response.data, indent=4))

        instance = DrugExNet.objects.get(pk=response.data["id"])
        builder_class = getattr(builders, builders.DrugExBuilder.__name__)
        builder = builder_class(
            instance
        )
        builder.build()



