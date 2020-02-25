import json

from django.test import TestCase
from django.urls import reverse
from rest_framework.test import APITestCase

from generators.apps import GeneratorsConfig
from generators.models import Generator, GeneratedMolSet
from modelling.apps import ModellingConfig
from projects.models import Project
from qsar.apps import QsarConfig


class MapTestCase(APITestCase):

    def createProject(self):
        post_data = {
          "name": "Test Project to POST",
          "description": "test description",
        }
        create_url = reverse('project-list')
        response = self.client.post(create_url, data=post_data, format='json')
        self.assertEqual(response.status_code, 201)
        print(json.dumps(response.data, indent=4))

        return Project.objects.get(pk=response.data["id"])

    def setUp(self):
        GeneratorsConfig.ready('dummy', True)
        ModellingConfig.ready('dummy', True)
        QsarConfig.ready('dummy', True)
        self.project = self.createProject()
        self.generator = Generator.objects.filter(project=self.project).all()[0]

        post_data = {
            "source" : self.generator.id,
            "name" : "Some Generated Set",
            "project" : self.project.id,
            "nSamples" : 100,
        }
        response = self.client.post(reverse('generatedSet-list'), data=post_data, format='json')
        self.assertEqual(response.status_code, 201)
        self.molset = GeneratedMolSet.objects.get(pk=response.data[1]["id"])

    def tearDown(self) -> None:
        self.project.delete()

    def test_create_map(self):
        # TODO: test this
        pass

