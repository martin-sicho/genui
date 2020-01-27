from rest_framework.test import APITestCase
from qsar.tests import setUp, createTestModel

class DrugExGeneratorInitTestCase(APITestCase):

    def setUp(self):
        from generators.apps import GeneratorsConfig
        GeneratorsConfig.ready('dummy')
        setUp(self)
        self.environ = createTestModel(self)

    def test_create_view(self):
        pass


