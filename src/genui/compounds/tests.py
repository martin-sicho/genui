"""
tests

Created by: Martin Sicho
On: 4/28/20, 2:05 PM
"""
import json

from genui.projects.tests import ProjectMixIn
from . import models

class CompoundsMixIn(ProjectMixIn):

    def createMolSet(self, url, appendData):
        post_data = {
            "name": "Test ChEMBL Data Set",
            "description": "Some description...",
            "project": self.project.id,

        }
        post_data.update(appendData)
        response = self.client.post(url, post_data)
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 201)

        return models.MolSet.objects.get(pk=response.data['id'])

