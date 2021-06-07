"""
tests

Created by: Martin Sicho
On: 4/28/20, 2:05 PM
"""
import json

from rest_framework.test import APITestCase

from genui.projects.tests import ProjectMixIn
from . import models
from .initializers.base import MolSetInitializer

class CompoundsMixIn(ProjectMixIn):

    def createMolSet(self, url, appendData):
        post_data = {
            "name": "Test Compounds",
            "description": "Some description...",
            "project": self.project.id,

        }
        post_data.update(appendData)
        response = self.client.post(url, post_data)
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 201)

        return models.MolSet.objects.get(pk=response.data['id'])

class DuplicateEntitiesTest(ProjectMixIn, APITestCase):

    class InitializerWithMultipleEntities(MolSetInitializer):

        def populateInstance(self) -> int:
            mol = self.addMoleculeFromSMILES('CCO')
            return mol.id

        def updateInstance(self) -> int:
            pass

    def test_create(self):
        molset = models.MolSet.objects.create(project=self.project, name='Test Set')
        entity = models.ChemicalEntity.objects.create(
            canonicalSMILES="CCO",
            inchi="1S/C2H6O/c1-2-3/h3H,2H2,1H3",
            inchiKey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        )

        mol_orig = models.Molecule(
            entity=entity
        )
        mol_orig.save()
        mol_orig.providers.add(molset)
        mol_orig.save()

        mol_duplicate = models.Molecule(
            entity=entity
        )
        mol_duplicate.save()
        mol_duplicate.providers.add(molset)
        mol_duplicate.save()

        initializer = self.InitializerWithMultipleEntities(molset)
        initializer.populateInstance()
