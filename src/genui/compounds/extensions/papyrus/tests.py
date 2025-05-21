import json
import os

from django.urls import reverse
from rest_framework.test import APITestCase
from genui.compounds.tests import CompoundsMixIn

from . import models

class PapyrusMolSetTestCase(CompoundsMixIn, APITestCase):

    def setUp(self):
        super().setUp()
        self.project = self.createProject()
        self.molset = self.createMolSet(
            reverse('papyrusSet-list'),
            {
                "targets": ["P29274", "P00533"], # P08908 - 5HT1A receptor
                "maxPerTarget" : 5
            }
        )
        self.assertEqual(self.molset.__class__, models.PapyrusCompounds)
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
    '''
    def test_no_max_per_target(self):
        self.molset = self.createMolSet(
            reverse('papyrusSet-list'),
            {
                "targets": ["O60885"],
            }
        )
    '''