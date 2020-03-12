import json

from django.test import TestCase

# Create your tests here.
from django.urls import reverse

from rest_framework.test import APITestCase

from generators.apps import GeneratorsConfig
from generators.models import Generator
from modelling.apps import ModellingConfig
from projects.models import Project


class ProjectMixIn:

    def createProject(self):
        post_data = {
          "name": "Test Project to POST",
          "description": "test description",
        }
        create_url = reverse('project-list')
        response = self.client.post(create_url, data=post_data, format='json')
        self.assertEqual(response.status_code, 201)
        print(json.dumps(response.data, indent=4))

        return Project.objects.get(pk=response.data['id'])

class ProjectTestCase(ProjectMixIn, APITestCase):

    def setUp(self) -> None:
        GeneratorsConfig.ready('dummy', True)
        ModellingConfig.ready('dummy', True)
        self.project = self.createProject()

    def test_default_generator(self):
        generator = Generator.objects.filter(project=self.project).all()[0]
        generator.get(100)

        generator.project.delete()
