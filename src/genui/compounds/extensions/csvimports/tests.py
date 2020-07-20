import json
import os

from django.urls import reverse
from rest_framework.test import APITestCase
from genui.projects.tests import ProjectMixIn


class CSVMolSetTestCase(ProjectMixIn, APITestCase):

    def setUp(self):
        super().setUp()
        self.project = self.createProject()
        self.filePath = os.path.join(os.path.dirname(__file__), 'test_files/init.csv')

    def test_csv_create(self):
        post_data = {
            'file': open(self.filePath),
            'name': "Molecule Set from a CSV Table",
            'project': self.project.id
        }
        url = reverse('csvSet-list')
        response = self.client.post(url, post_data)
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 201)

        # get the object detail from API
        url = reverse('csvSet-detail', args=[response.data['id']])
        response = self.client.get(url)
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 200)

        # get the activities
        url = reverse('activitySet-activities', args=[response.data['activities'][0]])
        response = self.client.get(url)
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 200)
        self.assertTrue(response.data['count'] == 3)
