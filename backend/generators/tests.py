import json
import os

from django.urls import reverse
from rest_framework.test import APITestCase

from generators.core import builders
from generators.models import DrugExNet
from modelling.models import Algorithm, AlgorithmMode
from qsar.tests import InitMixIn

class SetUpGeneratorsMixIn(InitMixIn):

    def createGenerator(self, create_url, initial=None):
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
        if initial:
            post_data["parent"] = initial
        response = self.client.post(create_url, data=post_data, format='json')
        self.assertEqual(response.status_code, 201)
        print(json.dumps(response.data, indent=4))

        instance = DrugExNet.objects.get(pk=response.data["id"])
        builder_class = getattr(builders, builders.DrugExBuilder.__name__)
        builder = builder_class(
            instance
        )
        builder.build()

        return instance

    def setUp(self):
        super().setUp()
        from generators.apps import GeneratorsConfig
        GeneratorsConfig.ready('dummy', True)
        self.drugex1 = self.createGenerator(reverse("drugex_net-list"))
        self.drugex2 = self.createGenerator(reverse("drugex_net-list"), initial=self.drugex1.id)


class DrugExGeneratorInitTestCase(SetUpGeneratorsMixIn, APITestCase):

    def setUp(self):
        super().setUp()
        self.environ = self.createTestModel()

    def test_create_drugexnet_view(self):
        self.assertTrue(self.drugex2.parent.id == self.drugex1.id)

        perf_view = reverse("drugex_perf_view", args=[self.drugex1.id])
        response = self.client.get(perf_view)
        self.assertEqual(response.status_code, 200)
        print(json.dumps(response.data, indent=4))
        self.assertTrue(response.data["count"] == 10)

        self.project.delete()
        self.assertFalse(os.path.exists(self.drugex1.modelFile.path))
        self.assertFalse(os.path.exists(self.drugex2.modelFile.path))


