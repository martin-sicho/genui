import json
import shutil

from django.conf import settings
from django.contrib.auth.models import User, Group
from django.core.management import call_command
from django.urls import reverse

from rest_framework.test import APITestCase

from genui.generators.models import Generator
from genui.projects.models import Project


class UserMixIn:

    def setUp(self):
        call_command('genuisetup', force=1, strict=0)
        self.user, self.username, self.password, *rest = self.createUser()
        self.login()

    def tearDown(self):
        # self.deleteUser()
        shutil.rmtree(settings.MEDIA_ROOT)

    def login(self):
        self.client.login(username=self.username, password=self.password)
        self.clientLoggedIn = True
        # self.clientSession = self.client.session
        # self.clientSession['CELERY_TASK_ALWAYS_EAGER'] = [settings.CELERY_TASK_ALWAYS_EAGER]
        # self.clientSession.save()

    def deleteUser(self):
        self.user.delete()
        self.clientLoggedIn = False
        self.username = None
        self.password = None


    def createUser(self):
        password = "temptestuser123456789"
        username = "temptestuser"
        email = "temptestuser@example.com"
        user = User.objects.create_user(username=username, email=email, password=password)
        usergroup = Group.objects.get(name='GenUI_Users')
        user.groups.add(usergroup)
        user.save()
        return user, username, password, email

class ProjectMixIn(UserMixIn):

    def setUp(self) -> None:
        super().setUp()
        self.project = self.createProject()

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

    def test_default_generator(self):
        self.assertTrue(self.project)
        generator = Generator.objects.filter(project=self.project).all()[0]
        self.assertTrue(generator)
