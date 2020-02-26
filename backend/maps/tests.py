import json

from django.test import TestCase
from django.urls import reverse
from rest_framework.test import APITestCase

from compounds.initializers import ChEMBLSetInitializer
from compounds.models import ChEMBLCompounds, Molecule
from generators.apps import GeneratorsConfig
from generators.models import Generator, GeneratedMolSet
from maps.apps import MapsConfig
from maps.core import builders
from maps.models import Map
from modelling.apps import ModellingConfig
from modelling.models import Algorithm, AlgorithmMode
from projects.models import Project
from qsar.apps import QsarConfig
from qsar.models import DescriptorGroup


class MapTestCase(APITestCase):

    def createProject(self):
        post_data = {
          "name": "Test Project to POST",
          "description": "test description",
        }
        create_url = reverse('project-list')
        response = self.client.post(create_url, data=post_data, format='json')
        self.assertEqual(response.status_code, 201)
        print(json.dumps(response.data, indent=4))

        return Project.objects.get(pk=response.data["id"])

    def createMolset(self, target):
        ret = ChEMBLCompounds.objects.create(**{
            "name": f"{target} Test",
            "description": "Some description...",
            "project": self.project
        })
        initializer = ChEMBLSetInitializer(
            ret
            , targets=[target]
            , max_per_target=10
        )
        initializer.populateInstance()
        return ret

    def setUp(self):
        ModellingConfig.ready('dummy', True)
        QsarConfig.ready('dummy', True)
        GeneratorsConfig.ready('dummy', True)
        MapsConfig.ready('dummy', True)
        self.project = self.createProject()
        self.molsets = [
            self.createMolset("CHEMBL251")
            , self.createMolset("CHEMBL203")
        ]

    def tearDown(self) -> None:
        self.project.delete()

    def test_create_map(self):
        post_data = {
            "name": "Test TSNE Map",
            "project": self.project.id,
            "trainingStrategy": {
                "algorithm": Algorithm.objects.get(name="TSNE").id,
                "mode": AlgorithmMode.objects.get(name="map").id,
                "parameters": {
                    "n_iter": 500,
                },
                "descriptors": [
                    DescriptorGroup.objects.get(name="MORGANFP").id
                ]
            },
            "molsets" : [x.id for x in self.molsets]
        }
        create_url = reverse('map-list')
        response = self.client.post(create_url, data=post_data, format='json')
        self.assertEqual(response.status_code, 201)
        print(json.dumps(response.data, indent=4))

        mymap = Map.objects.get(pk=response.data["id"])
        builder_class = getattr(builders, mymap.builder.name)
        builder = builder_class(
            mymap
        )
        builder.build()

        points_url = reverse('map-points-list', args=[mymap.id])
        response = self.client.get(points_url)
        self.assertEqual(response.status_code, 200)
        print(json.dumps(response.data, indent=4))
