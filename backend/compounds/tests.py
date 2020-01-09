from django.test import TestCase
from django.urls import reverse

from compounds.initializers.chembl import ChEMBLSetInitializer
from compounds.models import ChEMBLCompounds, MolSet
from projects.models import Project

class ChEMBLMolSetTestCase(TestCase):

    def setUp(self):
        self.project = Project.objects.create(**{
            "name" : "Test Project"
            , "description" : "Test Description"
        })

    def test_populate_task(self):
        molset = ChEMBLCompounds.objects.create(**{
            "name": "Test ChEMBL Data Set",
            "description": "Some description...",
            "project": self.project
        })
        instance = MolSet.objects.get(pk=molset.id)
        initializer = ChEMBLSetInitializer(instance, targets=["CHEMBL251", "CHEMBL203"], max_per_target=200)
        initializer.populateInstance()
        self.assertGreater(instance.molecules.count(), 0)
        for err in initializer.errors:
            self.assertTrue("Missing canonical SMILES string" in repr(err))
        molset.delete()


    # as of now this test affects the production database
    # TODO: create a separate database and celery worker to test this properly
    # def test_populate_view(self):
    #     populate_url = reverse('chemblSet-list')
    #     response = self.client.post(populate_url, {
    #         "name": "Test ChEMBL Data Set",
    #         "description": "Some description...",
    #         "project": self.project.id,
    #         "targets": [
    #             "CHEMBL251",
    #         ],
    #         "maxPerTarget" : 1000
    #     })
    #     self.assertEqual(response.status_code, 201)
    #     self.assertEqual(MolSet.objects.count(), 1)
    #     self.assertEqual(ChEMBLCompounds.objects.count(), 1)
