import json

from django.test import TestCase
from django.urls import reverse

from projects.tests import ProjectMixIn
from . import models, initializers

class CompoundsMixIn(ProjectMixIn):

    def getMolSet(self, targets, max_per_target=50):
        molset = models.ChEMBLCompounds.objects.create(**{
            "name": "Test ChEMBL Data Set",
            "description": "Some description...",
            "project": self.project
        })
        instance = models.MolSet.objects.get(pk=molset.id)
        initializer = initializers.chembl.ChEMBLSetInitializer(
            instance
            , targets=targets
            , max_per_target=max_per_target
        )
        initializer.populateInstance()
        self.assertGreater(instance.molecules.count(), 0)
        for err in initializer.errors:
            print(err)
            self.assertTrue("Missing canonical SMILES string" in repr(err) or "No standard activity value found for molecule" in repr(err))

        return molset

class ChEMBLMolSetTestCase(CompoundsMixIn, TestCase):

    def setUp(self):
        super().setUp()
        self.project = self.createProject()
        self.molset = self.getMolSet(["CHEMBL251", "CHEMBL203"], max_per_target=10)

    def test_activity_summary(self):
        activity_set = self.molset.activities.all()[0]
        summary_url = reverse('activitySet-summary', args=[activity_set.id])
        response = self.client.get(summary_url)
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 200)

    # def test_populate_view(self):
    #     populate_url = reverse('chemblSet-list')
    #     response = self.client.post(populate_url, {
    #         "name": "Test ChEMBL Data Set",
    #         "description": "Some description...",
    #         "project": self.project.id,
    #         "targets": [
    #             "CHEMBL251",
    #         ],
    #         "maxPerTarget" : 50
    #     })
    #     self.assertEqual(response.status_code, 201)
