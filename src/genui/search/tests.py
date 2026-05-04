import json

from rest_framework.test import APITestCase
from django.urls import reverse

from genui.compounds.extensions.chembl.tests import CompoundsMixIn
from genui.compounds.models import Molecule


URLS = {
    # MolSet-scoped search endpoints
    "molset_similarity": "molsets_sim_search",
    "molset_substructure": "molsets_sub_search",
    "molset_smarts": "molsets_smarts_search",

    # Project-scoped search endpoints
    "project_similarity": "projects_sim_search",
    "project_substructure": "projects_sub_search",
    "project_smarts": "projects_smarts_search",

    # InChIKey endpoints
    "inchikey_occurrence": "occurrence_inchikey_search",
    "inchikey_projects": "projects_inchikey_search",
}


class SearchMixIn(CompoundsMixIn):

    def setUp(self):
        super().setUp()
        self.project = self.createProject()
        self.molset = self.createMolSet(
            reverse("chemblSet-list"),
            {
                "targets": ["CHEMBL251"],
                "maxPerTarget": 10,
            },
        )

    def get_any_molecule_from_molset(self):
        """
        Tries to fetch a molecule that belongs to the created MolSet.
        Adjust if your relation differs.
        """
        mol = self.molset.molecules.select_related("entity").first()
        self.assertIsNotNone(mol, "MolSet has no molecules; cannot run search tests.")
        return mol

    def get_any_smiles(self):
        """
        Uses the stored canonical SMILES from your entity model if present.
        Falls back to something valid if not.
        """
        mol = self.get_any_molecule_from_molset()
        entity = getattr(mol, "entity", None)
        self.assertTrue(
            entity is not None and hasattr(entity, "canonicalSMILES") and entity.canonicalSMILES,
            "No SMILES found on molecule entity; cannot run similarity and substructure tests.",
        )
        return getattr(entity, "canonicalSMILES")


    def get_any_inchikey(self):
        mol = self.get_any_molecule_from_molset()
        entity = getattr(mol, "entity", None)
        self.assertTrue(
            entity is not None and hasattr(entity, "inchiKey") and entity.inchiKey,
            "No InChIKey found on molecule entity; cannot run inchikey tests.",
        )
        return entity.inchiKey
    

class SearchEndpointsTestCase(SearchMixIn, APITestCase):

    # Molset-based searches
    def test_molset_similarity_search(self):
        url = reverse(URLS["molset_similarity"])

        post_data = {
            "ids": [self.molset.id],
            "input": "c1ccccc1",
            "fp_type": "morganFP",
            "metric": "tanimoto",
            "threshold": 0,
            "top_n": 10,
        }

        response = self.client.post(url, data=post_data, format="json")

        self.assertEqual(response.status_code, 200)
        self.assertIn("query", response.data)
        self.assertIn("hits", response.data)
        self.assertIn("total_searched", response.data)
        self.assertIn("total_returned", response.data)

        self.assertEqual(response.data["query"]["ids"], [self.molset.id])
        self.assertTrue(isinstance(response.data["hits"], list))
        self.assertTrue(len(response.data["hits"]) > 0)
        self.assertGreaterEqual(response.data["total_searched"],response.data["total_returned"])

    def test_molset_substructure_search(self):
        url = reverse(URLS["molset_substructure"])

        post_data = {
            "ids": [self.molset.id],
            "input": "c1ccccc1",
        }

        response = self.client.post(url, data=post_data, format="json")

        self.assertEqual(response.status_code, 200)
        self.assertIn("query", response.data)
        self.assertIn("hits", response.data)
        self.assertIn("total_searched", response.data)
        self.assertIn("total_returned", response.data)

        self.assertEqual(response.data["query"]["ids"], [self.molset.id])
        self.assertTrue(isinstance(response.data["hits"], list))
        self.assertTrue(len(response.data["hits"]) > 0)
        self.assertGreaterEqual(response.data["total_searched"],response.data["total_returned"])


    def test_molset_smarts_search(self):
        url = reverse(URLS["molset_smarts"])

        post_data = {
            "ids": [self.molset.id],
            "input": "c1ccccc1",
        }

        response = self.client.post(url, data=post_data, format="json")

        self.assertEqual(response.status_code, 200)
        self.assertIn("query", response.data)
        self.assertIn("hits", response.data)
        self.assertIn("total_searched", response.data)
        self.assertIn("total_returned", response.data)

        self.assertEqual(response.data["query"]["ids"], [self.molset.id])
        self.assertTrue(isinstance(response.data["hits"], list))
        self.assertTrue(len(response.data["hits"]) > 0)
        self.assertGreaterEqual(response.data["total_searched"],response.data["total_returned"])
    
    def test_molset_search_returns_404_for_missing_molset(self):
        url = reverse(URLS["molset_similarity"])

        post_data = {
            "ids": [99999999999],
            "input": "c1ccccc1",
            "fp_type": "morganFP",
            "metric": "tanimoto",
            "threshold": 0,
            "top_n": 5,
        }

        response = self.client.post(url, data=post_data, format="json")

        self.assertEqual(response.status_code, 404)
        self.assertIn("error", response.data)

    def test_number_of_hits(self):
        url = reverse(URLS["molset_similarity"])
        n_compounds = len(self.molset.molecules.select_related("entity").all())

        post_data = {
            "ids": [self.molset.id],
            "input": "c1ccccc1",
            "fp_type": "morganFP",
            "metric": "tanimoto",
            "threshold": 0,
            "top_n": n_compounds + 5,
        }

        response = self.client.post(url, data=post_data, format="json")

        self.assertTrue(len(response.data["hits"]) == min(post_data["top_n"], n_compounds))

    def test_unique_compounds(self):
        url = reverse(URLS["molset_similarity"])

        post_data = {
            "ids": [self.molset.id],
            "input": self.get_any_smiles(),
            "fp_type": "morganFP",
            "metric": "tanimoto",
            "threshold": 1,
            "top_n": 5,
        }

        response = self.client.post(url, data=post_data, format="json")
        print(json.dumps(response.data, indent=4))
        self.assertTrue(len(response.data["hits"]) == 1)

    # Project-based searches
    def test_project_similarity_search(self):
        url = reverse(URLS["project_similarity"])

        post_data = {
            "ids": [self.project.id],
            "input": "c1ccccc1",
            "fp_type": "morganFP",
            "metric": "tanimoto",
            "threshold": 0,
            "top_n": 10,
        }

        response = self.client.post(url, data=post_data, format="json")

        self.assertEqual(response.status_code, 200)
        self.assertIn("query", response.data)
        self.assertIn("hits", response.data)
        self.assertIn("total_searched", response.data)
        self.assertIn("total_returned", response.data)

        self.assertEqual(response.data["query"]["ids"], [self.project.id])
        self.assertTrue(isinstance(response.data["hits"], list))
        self.assertTrue(len(response.data["hits"]) > 0)
        self.assertGreaterEqual(response.data["total_searched"],response.data["total_returned"])

    def test_project_substructure_search(self):
        url = reverse(URLS["project_substructure"])

        post_data = {
            "ids": [self.project.id],
            "input": "c1ccccc1",
        }

        response = self.client.post(url, data=post_data, format="json")

        self.assertEqual(response.status_code, 200)
        self.assertIn("query", response.data)
        self.assertIn("hits", response.data)
        self.assertIn("total_searched", response.data)
        self.assertIn("total_returned", response.data)

        self.assertEqual(response.data["query"]["ids"], [self.project.id])
        self.assertTrue(isinstance(response.data["hits"], list))
        self.assertTrue(len(response.data["hits"]) > 0)
        self.assertGreaterEqual(response.data["total_searched"],response.data["total_returned"])

    def test_project_smarts_search(self):
        url = reverse(URLS["project_smarts"])

        post_data = {
            "ids": [self.project.id],
            "input": "c1ccccc1[O,S]",
        }

        response = self.client.post(url, data=post_data, format="json")

        self.assertEqual(response.status_code, 200)
        self.assertIn("query", response.data)
        self.assertIn("hits", response.data)
        self.assertIn("total_searched", response.data)
        self.assertIn("total_returned", response.data)

        self.assertEqual(response.data["query"]["ids"], [self.project.id])
        self.assertTrue(isinstance(response.data["hits"], list))
        self.assertTrue(len(response.data["hits"]) > 0)
        self.assertGreaterEqual(response.data["total_searched"],response.data["total_returned"])

    def test_project_search_returns_404_for_missing_project(self):
        url = reverse(URLS["project_similarity"])

        post_data = {
            "ids": [999999999],
            "input": "c1ccccc1",
            "fp_type": "morganFP",
            "metric": "tanimoto",
            "threshold": 0,
            "top_n": 5,
        }

        response = self.client.post(url, data=post_data, format="json")

        self.assertEqual(response.status_code, 404)
        self.assertIn("error", response.data)


    # InChIKey endpoints
    def test_inchikey_projects_search(self):
        url = reverse(URLS["inchikey_projects"])

        post_data = {
            "input": self.get_any_inchikey(),
        }

        response = self.client.post(url, data=post_data, format="json")

        self.assertEqual(response.status_code, 200)
        self.assertIn("query", response.data)
        self.assertIn("hits", response.data)
        self.assertIn("total_searched", response.data)
        self.assertIn("total_returned", response.data)

        self.assertTrue(isinstance(response.data["hits"], list))
        self.assertTrue(len(response.data["hits"]) == 1)
        self.assertGreaterEqual(response.data["total_searched"],response.data["total_returned"])

    def test_inchikey_occurrence_search(self):
        url = reverse(URLS["inchikey_occurrence"])

        post_data = {
            "input": self.get_any_inchikey(),
        }

        response = self.client.post(url, data=post_data, format="json")

        self.assertEqual(response.status_code, 200)
        self.assertIn("query", response.data)
        self.assertIn("occurrence", response.data)