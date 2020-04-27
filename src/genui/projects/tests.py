import json

from django.contrib.auth.models import User

# Create your tests here.
from django.urls import reverse

from rest_framework.test import APITestCase

from genui.compounds.apps import CompoundsConfig
from genui.generators.apps import GeneratorsConfig
from genui.generators.models import Generator
from genui.maps.apps import MapsConfig
from genui.modelling.apps import ModellingConfig
from genui.projects.apps import ProjectsConfig
from genui.projects.models import Project
from genui.qsar.apps import QsarConfig


class UserMixIn:

    def setUp(self):
        self.user, self.username, self.password, *rest = self.createUser()
        self.login()

    def login(self):
        self.client.login(username=self.username, password=self.password)

    def createUser(self):
        password = "temptestuser123456789"
        username = "temptestuser"
        email = "temptestuser@example.com"
        return User.objects.create_user(username=username, email=email, password=password), username, password, email

class ProjectMixIn(UserMixIn):

    def setUp(self):
        super().setUp()
        ProjectsConfig.ready('dummy', True)
        CompoundsConfig.ready('dummy', True)
        ModellingConfig.ready('dummy', True)
        GeneratorsConfig.ready('dummy', True)
        QsarConfig.ready('dummy', True)
        MapsConfig.ready('dummy', True)

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
        super().setUp()
        self.project = self.createProject()

    def test_default_generator(self):
        generator = Generator.objects.filter(project=self.project).all()[0]
        print(generator.get(100))

        generator.project.delete()
