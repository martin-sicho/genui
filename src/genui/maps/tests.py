import json

from django.urls import reverse
from rest_framework.test import APITestCase

from genui.compounds.tests import CompoundsMixIn
from genui.maps.core import builders
from genui.maps.models import Map
from genui.modelling.models import Algorithm, AlgorithmMode
from genui.qsar.models import DescriptorGroup


class MapTestCase(CompoundsMixIn, APITestCase):

    def setUp(self):
        super().setUp()
        self.project = self.createProject()
        self.molsets = [
            self.getMolSet(["CHEMBL251"], max_per_target=10),
            self.getMolSet(["CHEMBL203"], max_per_target=10),
        ]

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
