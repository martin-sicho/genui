import json
import os

from django.urls import reverse
from rest_framework.test import APITestCase
from genui.compounds.tests import CompoundsMixIn

from . import models

class ChEMBLMolSetTestCase(CompoundsMixIn, APITestCase):

    def setUp(self):
        super().setUp()
        self.project = self.createProject()
        self.molset = self.createMolSet(
            reverse('chemblSet-list'),
            {
                "targets": ["CHEMBL251", "CHEMBL203"],
                "maxPerTarget" : 10
            }
        )
        self.assertEquals(self.molset.__class__, models.ChEMBLCompounds)
        mol = self.molset.molecules.all()[0]
        image_path = mol.mainPic.image.path
        self.assertTrue(os.path.exists(image_path))
        mol.delete()
        self.assertFalse(os.path.exists(image_path))

    def test_activity_summary(self):
        activity_set = self.molset.activities.all()[0]
        summary_url = reverse('activitySet-summary', args=[activity_set.id])
        response = self.client.get(summary_url)
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 200)
