import json

from django.urls import reverse
from rest_framework.test import APITestCase

from genui.compounds.models import MolSet
from genui.generators.extensions.genuidrugex.genuimodels import builders
from genui.compounds.extensions.generated.models import GeneratedMolSet
from genui.generators.extensions.genuidrugex.models import DrugExNet, DrugExAgent, ScoringMethod, ScoreModifier, \
    DrugExEnvironment, DrugExScorer, DrugExAgentTraining, DrugEx, DrugExNetTraining
from genui.models.models import Algorithm, AlgorithmMode, ModelFile, Model
from genui.qsar.tests import QSARModelInit

TEST_EPOCHS = 2

class SetUpDrugExGeneratorsMixIn(QSARModelInit):

    def getPerformance(self, url):
        response = self.client.get(url)
        self.assertEqual(response.status_code, 200)
        print(json.dumps(response.data, indent=4))
        return response

    def createDrugExNet(self, create_url, initial=None):
        post_data = {
          "name": "Test DrugEx Network (pretraining)" if not initial else "Test DrugEx Network (finetuning)",
          "description": "test description",
          "project": self.project.id,
          "trainingStrategy": {
            "algorithm": Algorithm.objects.get(name="DrugExNetwork").id,
            "mode": AlgorithmMode.objects.get(name="generator").id,
            "parameters": {
                "nEpochs": TEST_EPOCHS,
                "batchSize" : 16,
            },
            "modelClass": 'SS',
            "inputType" : 'MS'
          },
          "validationStrategy": {
            "validSetSize": 5
          },
          "molset": self.molset.id
        }
        if initial:
            post_data["parent"] = initial.id
        response = self.client.post(create_url, data=post_data, format='json')
        self.assertEqual(response.status_code, 201)
        print(json.dumps(response.data, indent=4))

        return DrugExNet.objects.get(pk=response.data["id"])

    def createDrugExAgent(self, url, exploit_net, explore_net, environ):
        post_data = {
            "name": "Test DrugEx Agent",
            "description": "test description",
            "project": self.project.id,
            "trainingStrategy": {
                "algorithm": Algorithm.objects.get(name="DrugExAgent").id,
                "mode": AlgorithmMode.objects.get(name="generator").id,
                "explorer": DrugExAgentTraining.ExplorerClass.smiles_molecules,
                "parameters": {
                    "nEpochs": TEST_EPOCHS,
                    "batchSize" : 16,
                    "epsilon" : 0.01,
                    "beta" : 0.1,
                },
            },
            "validationStrategy": {
                "validSetSize": 5
            },
            "environment": environ.id,
            "exploitationNet" : exploit_net.id,
            "explorationNet" : explore_net.id,
        }
        response = self.client.post(url, data=post_data, format='json')
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 201)

        return DrugExAgent.objects.get(pk=response.data["id"])

    def createModelScorer(self, model):
        response = self.client.post(
            reverse('drugex_scoremethods_genuimodels-list'),
            data = {
                "name": "Test GenUI QSAR Model Scorer",
                "description": "scorer test description",
                "project": self.project.id,
                "model": model.id
            },
            format='json'
        )
        self.assertEqual(response.status_code, 201)
        print(json.dumps(response.data, indent=4))

        return ScoringMethod.objects.get(pk=response.data["id"])

    def createPropertyScorer(self, prop):
        response = self.client.post(
            reverse('drugex_scoremethods_properties-list'),
            data = {
                "name": f"Test GenUI {prop} Property Scorer",
                "description": "scorer test description",
                "project": self.project.id,
                "prop": prop
            },
            format='json'
        )
        self.assertEqual(response.status_code, 201)
        print(json.dumps(response.data, indent=4))

        return ScoringMethod.objects.get(pk=response.data["id"])

    def createScoreModifier(self, url, settings):
        data = {
            "name": f"Test GenUI Score Modifier",
            "description": "modifier test description",
            "project": self.project.id,
            **settings
        }
        response = self.client.post(
            url,
            data = data,
            format='json'
        )
        self.assertEqual(response.status_code, 201)
        print(json.dumps(response.data, indent=4))

        return ScoreModifier.objects.get(pk=response.data["id"])

    def createScorer(self, environment, method, modifier, threshold):
        url = reverse('drugex_scorers-list')
        response = self.client.post(
            url,
            data = {
                'name': 'Test Scorer for Drugex',
                'description': 'scorer test description',
                'project': self.project.id,
                'environment': environment.id,
                'modifier': modifier.id,
                'method': method.id,
                'threshold': threshold
            },
            format='json'
        )
        self.assertEqual(response.status_code, 201)
        print(json.dumps(response.data, indent=4))

        return DrugExScorer.objects.get(pk=response.data['id'])

    def createTestEnvironment(self):
        # scoring functions
        model_scorer = self.createModelScorer(
            self.createTestQSARModel()
        )
        property_scorer_sa = self.createPropertyScorer("SA")
        property_scorer_qed = self.createPropertyScorer("QED")

        # modifiers
        clipped_modifier_url = reverse('drugex_modifiers_clipped-list')
        modifier_qed = self.createScoreModifier(
            clipped_modifier_url,
            {
                'lower': 0,
                'upper': 1
            }
        )
        modifier_sa = self.createScoreModifier(
            clipped_modifier_url,
            {
                'lower': 4.5,
                'upper': 0
            }
        )
        modifier_qsar = self.createScoreModifier(
            clipped_modifier_url,
            {
                'lower': 0.2,
                'upper': 0.5
            }
        )

        # create the environment
        response = self.client.post(
            reverse('drugex_env-list'),
            data={
                'project' : self.project.id,
                'name' : 'Test DrugEx Environment',
                'description': 'environment description for testing',
                'rewardScheme': DrugExEnvironment.RewardScheme.paretoCrowding
            },
            format='json'
        )
        self.assertEqual(response.status_code, 201)
        print(json.dumps(response.data, indent=4))
        environment = DrugExEnvironment.objects.get(pk=response.data['id'])

        # add scorers
        self.createScorer(
            environment,
            model_scorer,
            modifier_qsar,
            threshold=0.5
        )
        self.createScorer(
            environment,
            property_scorer_qed,
            modifier_qed,
            threshold=0.0
        )
        self.createScorer(
            environment,
            property_scorer_sa,
            modifier_sa,
            threshold=0.0
        )

        return environment

    def postDrugExModelFile(self, instance, data):
        url = reverse('drugex-net-model-files-list', args=[instance.id])
        response = self.client.post(
            url,
            data=data,
            format='multipart'
        )
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 201)

    def getBuilder(self, instance : Model, builder_module):
        builder_class = getattr(builder_module, instance.builder.name)
        return builder_class(
            instance
        )

    def createGenerator(self, agent, molset):
        generator_data = {
            "name": f"Test Generator for {agent.name}",
            "project": self.project.id,
            "molset": self.molset.id,
            "agent": agent.id
        }
        response = self.client.post(
            reverse('drugex_generators-list'),
            data=generator_data,
            format='json'
        )
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 201)
        generator = DrugEx.objects.get(pk=response.data['id'])

        response = self.client.get(reverse('generator-list'))
        self.assertEqual(response.status_code, 200)
        generator_found = False
        for item in response.data:
            if item['id'] == generator.id:
                generator_found = True
                mols = generator.get(100)
                self.assertTrue(len(mols) > 0 and type(mols[0]) == str)
        self.assertTrue(generator_found)

        return generator

    def createGeneratedMolSet(self, generator):
        generated_set_data = {
            "source" : generator.id,
            "name" : f"Test Generated DrugEx Set with {generator.name}",
            "project" : self.project.id,
            "nSamples" : 100,
        }
        response = self.client.post(reverse('generatedSet-list'), data=generated_set_data, format='json')
        self.assertEqual(response.status_code, 201)
        print(json.dumps(response.data, indent=4))
        mols = MolSet.objects.get(pk=response.data['id'])
        smiles = mols.allSmiles
        self.assertTrue(len(smiles) > 0 and type(smiles[0]) == str)

        return mols

class DrugExFromFileTestCase(SetUpDrugExGeneratorsMixIn, APITestCase):

    def test_create_from_files(self):
        instance_first = self.createDrugExNet(reverse("drugex_net-list"))
        upload_main = open(instance_first.modelFile.path, "rb")
        upload_voc = open(instance_first.files.filter(note=DrugExNet.VOC_FILE_NOTE).get().path, "rb")

        create_url = reverse('drugex_net-list')
        post_data = {
            "name": "Test DrugEx Net from File",
            "project": self.project.id,
            "build" : False,
            "trainingStrategy": {
                "algorithm": Algorithm.objects.get(name="DrugExNetwork").id,
                "mode": AlgorithmMode.objects.get(name="generator").id,
                "modelClass": DrugExNetTraining.ModelClass.graphTrans,
                "inputType":  DrugExNetTraining.StructureInputType.frags
            },
        }
        response = self.client.post(create_url, data=post_data, format='json')
        print(json.dumps(response.data, indent=4))
        self.assertEqual(response.status_code, 201)
        instance = DrugExNet.objects.get(pk=response.data["id"])
        self.assertFalse(instance.modelFile)

        self.postDrugExModelFile(instance, {
                "file" : upload_main,
                "kind": ModelFile.MAIN,
        })
        self.postDrugExModelFile(instance, {
            "file" : upload_voc,
            "note" : DrugExNet.VOC_FILE_NOTE
        })

        generator = self.createGenerator(instance, self.molset)

        self.createGeneratedMolSet(generator)

        # test correct exception raised when a build attempt is made
        builder = self.getBuilder(instance, builders)
        self.assertRaises(NotImplementedError, builder.build)
        # dataset = instance.getDataSetFromMolset(generator.molset, tempfile.NamedTemporaryFile().name)
        # model = builder.model
        # mols = model.sample(100, from_inputs=dataset)
        # self.assertTrue(len(mols) == 2)
        # self.assertTrue(len(mols[0]) == len(mols[1]))
        # self.assertTrue(len(mols[0]) > 0)
        # print(mols)

class DrugExGeneratorInitTestCase(SetUpDrugExGeneratorsMixIn, APITestCase):

    def test_all(self):
        drugex1 = self.createDrugExNet(reverse("drugex_net-list"))
        drugex2 = self.createDrugExNet(reverse("drugex_net-list"), initial=drugex1)
        self.assertTrue(drugex2.parent.id == drugex1.id)

        environ = self.createTestEnvironment()
        agent = self.createDrugExAgent(reverse("drugex_agent-list"), drugex1, drugex2, environ)


        response = self.getPerformance(reverse("drugex_net_perf_view", args=[drugex1.id]))
        self.assertTrue(response.data["count"] > 0)

        response = self.getPerformance(reverse("drugex_agent_perf_view", args=[agent.id]))
        self.assertTrue(response.data["count"] > 0)

        response = self.client.get(reverse('generator-list'))
        self.assertEqual(response.status_code, 200)
        print(json.dumps(response.data, indent=4))

        response = self.client.get(reverse('generator-list'))
        self.assertEqual(response.status_code, 200)
        print(json.dumps(response.data, indent=4))

        for generator in [drugex1, drugex2, agent]:
            generator = DrugEx.objects.get(agent=generator)
            post_data = {
                "source" : generator.id,
                "name" : f"Test Generated DrugEx Set with {generator.name}",
                "project" : self.project.id,
                "nSamples" : 100,
            }
            response = self.client.post(reverse('generatedSet-list'), data=post_data, format='json')
            print(json.dumps(response.data, indent=4))
            self.assertEqual(response.status_code, 201)
            molset_id = response.data['id']

            response = self.client.get(reverse('molset-list'))
            self.assertEqual(response.status_code, 200)
            print(json.dumps(response.data, indent=4))

            generated_set = [x for x in response.data if x['id'] == molset_id]
            self.assertTrue(len(generated_set) == 1)
            generated_set = GeneratedMolSet.objects.get(pk=generated_set[0]['id'])

            mols_url = reverse('moleculesInSet', args=[generated_set.id])
            response = self.client.get(mols_url)
            self.assertEqual(response.status_code, 200)
            print(json.dumps(response.data, indent=4))
            self.assertTrue(response.data['count'] >= 0)

class UseDefaultNetTestCase(SetUpDrugExGeneratorsMixIn, APITestCase):

    def test_all(self):
        models_2_test = [
            # self.project.model_set.get(name__contains='ZINC'),
            # self.project.model_set.get(name__contains='ChEMBL'),
            self.project.model_set.get(name__contains='Graph-Based Transformer for ChEMBL 27'),
        ]
        for parent in models_2_test:
            # finetuning
            explore_net = self.createDrugExNet(
                reverse("drugex_net-list"),
                initial=parent
            )

            # RL
            environ = self.createTestEnvironment()
            agent = self.createDrugExAgent(reverse("drugex_agent-list"), exploit_net=parent, explore_net=explore_net, environ=environ)

            # generating mols
            generator = self.createGenerator(agent, self.molset)
            mols = self.createGeneratedMolSet(generator)
            print(mols.allSmiles)

            # predict with environment
            url = reverse('drugex_env-calculate', args=(environ.id,))
            response = self.client.post(
                url,
                data = {
                    "molsets": [self.molset.id, mols.id],
                    "useModifiers": False
                },
                format='json'
            )
            self.assertEqual(response.status_code, 201)
            print(json.dumps(response.data, indent=4))

            for molset_id in response.data['molsets']:
                activities_ids = self.client.get(
                    reverse('molset-detail', args=[molset_id])
                ).data['activities']
                for activity_id in activities_ids:
                    activities = self.client.get(reverse('activitySet-activities', args=[activity_id])).data
                    self.assertTrue(activities["count"] > 0)
