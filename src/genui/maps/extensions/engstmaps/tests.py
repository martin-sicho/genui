from django.test import TestCase
import json

from django.urls import reverse
from rest_framework.test import APITestCase

from genui.compounds.extensions.chembl.tests import CompoundsMixIn
from genui.compounds.models import MolSet
from genui.maps.models import Map
from genui.models.models import Algorithm, AlgorithmMode
from genui.qsar.models import DescriptorGroup

class MapTestCase(CompoundsMixIn, APITestCase):

    def setUp(self) -> None:
        super().setUp()
        self.project = self.createProject()
        self.molsets = [
            self.createMolSet(
                reverse('chemblSet-list'),
                {
                    "targets": ["CHEMBL251"],
                    "maxPerTarget": 10
                }
            ),
            self.createMolSet(
                reverse('chemblSet-list'),
                {
                    "targets": ["CHEMBL203"],
                    "maxPerTarget": 10
                }
            ),
        ]


    def test_create_get_map_correct(self):
        post_data = {
            "name": "PCA test",
            "project": self.project.id,
            "trainingStrategy": {
                "algorithm": Algorithm.objects.get(name="PCA").id,
                "mode": AlgorithmMode.objects.get(name="map").id,
                "descriptors": [
                    DescriptorGroup.objects.get(name="MORGANFP").id
                ]
            },
            "molsets" : [x.id for x in self.molsets]
        }
        mol_count = 0
        for molset in MolSet.objects.filter(id__in=post_data["molsets"]).all():
            count = molset.molecules.count()
            mol_count += count

        create_url = reverse('map-list')
        response = self.client.post(create_url, data=post_data, format='json')
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 201)

        mymap = Map.objects.get(pk=response.data["id"])

        points_url = reverse('map-points-list', args=[mymap.id])
        print(points_url)
        response = self.client.get(points_url)
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 200)
        self.assertEqual(mol_count, len(response.data["results"]))

        maps_url = reverse('map-list')
        map_list_response = self.client.get(maps_url)
        print(json.dumps(map_list_response.data, indent=4))
        self.assertEqual(map_list_response.status_code, 200)
        self.assertEqual(map_list_response.data[0]["name"], "PCA test")



        map_id = map_list_response.data[0]["id"]
        maps_url = reverse('map-detail', args=[mymap.id])
        print(maps_url)
        map_response = self.client.get(maps_url)
        print(json.dumps(map_response.data, indent=4))
        self.assertEqual(map_response.status_code, 200)
        self.assertEqual(map_response.data["name"], "PCA test")




    def test_create_map_wrong(self):
        post_data = {
            "name": "PCA test wrong",
            "project": self.project.id,
            "molsets": [x.id for x in self.molsets]
        }

        create_url = reverse('map-list')
        response = self.client.post(create_url, data=post_data, format='json')
        print(json.dumps(response.data, indent=4))
        self.assertGreater(response.status_code, 399)

