import uuid

import numpy as np
from django.core.files.base import ContentFile
from django.db import models
from django.utils.translation import gettext_lazy as _
from drugex.data.corpus.corpus import SequenceCorpus
from drugex.data.corpus.vocabulary import VocGraph, VocGPT, VocSmiles
from drugex.data.datasets import GraphFragDataSet, SmilesFragDataSet, SmilesDataSet
from drugex.data.fragments import FragmentCorpusEncoder, GraphFragmentEncoder, FragmentPairsSplitter, \
    SequenceFragmentEncoder
from drugex.data.processing import Standardization, CorpusEncoder, RandomTrainTestSplitter
from drugex.molecules.converters.fragmenters import Fragmenter
from drugex.training import environment
from drugex.training.interfaces import Scorer
from drugex.training.models import GraphModel, GPT2Model, RNN
from drugex.training.models.explorer import GraphExplorer, SmilesExplorer
from drugex.training.rewards import ParetoCrowdingDistance, ParetoSimilarity, WeightedSum
from drugex.training.scorers import modifiers

from drugex.training.scorers.properties import Property

from genui.compounds.models import MolSet, ActivitySet
from genui.generators.extensions.genuidrugex.torchutils import cleanup
from genui.generators.models import Generator
from genui.models.models import Model, ModelFile, ValidationStrategy, TrainingStrategy, ModelPerfomanceNN
from genui.projects.models import DataSet
from genui.qsar.genuimodels.builders import BasicQSARModelBuilder
from genui.qsar.models import QSARModel
from genui.utils.models import OverwriteStorage


class DrugExNet(Model):
    CORPUS_FILE_TRAIN_NOTE = "drugex_corpus_train"
    CORPUS_FILE_VALID_NOTE = "drugex_corpus_test"
    VOC_FILE_NOTE = "drugex_voc"

    molset = models.ForeignKey(MolSet, on_delete=models.CASCADE, null=True)
    parent = models.ForeignKey("self", on_delete=models.CASCADE, null=True)

    def createCorpusFile(self, note, name):
        ret = self.files.filter(kind=ModelFile.AUXILIARY, note=note)
        if not ret:
            ret = ModelFile.create(
                self,
                name,
                ContentFile(''),
                note=note
            )
        else:
            ret = ret.get()
        return ret

    @property
    def corpusFileTrain(self):
        return self.createCorpusFile(DrugExNet.CORPUS_FILE_TRAIN_NOTE, "corpus_train.tsv")

    @property
    def corpusFileTest(self):
        return self.createCorpusFile(DrugExNet.CORPUS_FILE_VALID_NOTE, "corpus_test.tsv")

    def getDefaultVoc(self):
        modelClass = self.trainingStrategy.modelClass
        inputType = self.trainingStrategy.inputType
        if modelClass == "GT" and inputType == "FS":
            return VocGraph()
        elif modelClass == "ST" and inputType == "FS":
            return VocGPT(VocSmiles.defaultWords)
        elif modelClass == "SS" and inputType == "MS":
            return VocSmiles(VocSmiles.defaultWords)
        else:
            NotImplementedError(f"Vocabulary for model with these settings is unavailable: modelClass={modelClass} and inputType={inputType}")

    @property
    def vocFile(self):
        ret = self.files.filter(kind=ModelFile.AUXILIARY, note=self.VOC_FILE_NOTE)
        if not ret:
            ret = ModelFile.create(
                self,
                "voc.txt",
                ContentFile(''),
                note=DrugExNet.VOC_FILE_NOTE
            )

            if self.parent:
                vocabulary = self.parent.getDefaultVoc().fromFile(self.parent.vocFile.path)
            else:
                vocabulary = self.getDefaultVoc()
            vocabulary.toFile(ret.path)
        else:
            ret = ret.get()
        return ret

    def getCorpus(self, _file):
        modelClass = self.trainingStrategy.modelClass
        inputType = self.trainingStrategy.inputType
        if modelClass == "GT" and inputType == "FS":
            data = GraphFragDataSet(_file.path, rewrite=False)
            data.setVoc(VocGraph.fromFile(self.vocFile.path))
            return data
        elif modelClass == "ST" and inputType == "FS":
            data = SmilesFragDataSet(_file.path, rewrite=False)
            data.setVoc(VocGPT.fromFile(self.vocFile.path))
            return data
        elif modelClass == "SS" and inputType == "MS":
            data = SmilesDataSet(_file.path, rewrite=False)
            data.setVoc(VocSmiles.fromFile(self.vocFile.path))
            return data
        else:
            raise NotImplementedError(f"Corpus for this model configuration does not exist: modelClass={modelClass} and inputType={inputType}")

    @property
    def corpusTrain(self):
        return self.getCorpus(self.corpusFileTrain)

    @property
    def corpusTest(self):
        return self.getCorpus(self.corpusFileTest)

    @staticmethod
    def standardize(smiles, n_proc=4):
        standardizer = Standardization(n_proc=n_proc)
        return standardizer.apply(smiles)

    def getDataSetFromMolset(self, molset, data_file, rewrite=True, n_proc=4):
        modelClass = self.trainingStrategy.modelClass
        inputType = self.trainingStrategy.inputType

        if modelClass != "SS" and not molset:
            molset = self.molset
            if not molset:
                raise RuntimeError(f"Could not determine molecule set to create input for generative agent: {repr(self)}.")

        smiles = self.standardize(molset.allSmiles, n_proc=n_proc)
        if modelClass == "GT" and inputType == "FS":
            if rewrite:
                encoder = FragmentCorpusEncoder(
                    fragmenter=Fragmenter(4, 4, 'brics'),
                    encoder=GraphFragmentEncoder(
                        VocGraph.fromFile(self.vocFile.path, n_frags=4)
                    ),
                    n_proc=n_proc
                )
                out_data = GraphFragDataSet(data_file.path, rewrite=True)
                encoder.apply(smiles, encodingCollectors=[out_data])
                return out_data, molset
            else:
                return GraphFragDataSet(data_file.path, rewrite=False), molset
        elif modelClass == "ST" and inputType == "FS":
            if rewrite:
                encoder = FragmentCorpusEncoder(
                    fragmenter=Fragmenter(4, 4, 'brics'),
                    encoder=SequenceFragmentEncoder(
                        VocGPT.fromFile(self.vocFile.path)
                    ),
                    n_proc=n_proc
                )
                out_data = SmilesFragDataSet(data_file.path, rewrite=True)
                encoder.apply(smiles, encodingCollectors=[out_data])
                return out_data, molset
            else:
                return SmilesFragDataSet(data_file.path, rewrite=False), molset
        elif modelClass == "SS" and inputType == "MS":
            if rewrite:
                encoder = CorpusEncoder(
                    SequenceCorpus,
                    {
                      'vocabulary' : VocSmiles.fromFile(self.vocFile.path)
                    },
                    n_proc=n_proc
                )
                out_data = SmilesDataSet(data_file.path, rewrite=True)
                encoder.apply(smiles, out_data)
                return out_data, molset
            else:
                return SmilesDataSet(data_file.path, rewrite=False), molset
        else:
             raise NotImplementedError(f"Data set encoding for this model configuration is not available: modelClass={modelClass} and inputType={inputType}.")

    def prepareData(self, vocabulary=None, n_proc=4, chunk_size=1000):
        if self.molset:
            train = self.corpusTrain
            test = self.corpusTest
            vocabulary = vocabulary if vocabulary else train.getVoc() + test.getVoc()

            modelClass = self.trainingStrategy.modelClass
            inputType = self.trainingStrategy.inputType
            if modelClass == "GT" and inputType == "FS":
                encoder = FragmentCorpusEncoder(
                    fragmenter=Fragmenter(4, 4, 'brics'),
                    encoder=GraphFragmentEncoder(
                        vocabulary
                    ),
                    pairs_splitter=FragmentPairsSplitter(0.1, self.validationStrategy.validSetSize, make_unique=False),
                    n_proc=n_proc,
                    chunk_size=chunk_size
                )

                smiles = self.molset.allSmiles
                smiles = self.standardize(smiles, n_proc=n_proc)
                encoder.apply(smiles, encodingCollectors=[test, train])
                vocabulary = test.getVoc() + train.getVoc()
                vocabulary.toFile(self.vocFile.path)
            elif modelClass == "ST" and inputType == "FS":
                encoder = FragmentCorpusEncoder(
                    fragmenter=Fragmenter(4, 4, 'brics'),
                    encoder=SequenceFragmentEncoder(
                        vocabulary
                    ),
                    pairs_splitter=FragmentPairsSplitter(0.1, self.validationStrategy.validSetSize, make_unique=False),
                    n_proc=n_proc,
                    chunk_size=chunk_size
                )

                smiles = self.molset.allSmiles
                smiles = self.standardize(smiles, n_proc=n_proc)
                encoder.apply(smiles, encodingCollectors=[test, train])
                vocabulary = test.getVoc() + train.getVoc()
                vocabulary.toFile(self.vocFile.path)
            elif modelClass == "SS" and inputType == "MS":
                smiles = self.molset.allSmiles
                splitter = RandomTrainTestSplitter(0.1, self.validationStrategy.validSetSize)
                train_smiles, test_smiles = splitter(smiles)
                for smis, out_data in zip([train_smiles, test_smiles], [train, test]):
                    encoder = CorpusEncoder(
                        SequenceCorpus,
                        {
                          'vocabulary' : vocabulary
                        },
                        n_proc=n_proc
                    )
                    encoder.apply(smis, out_data)
                vocabulary.toFile(self.vocFile.path)
            else:
                raise NotImplementedError(f"Data set encoding for this model configuration is not available: modelClass={modelClass} and inputType={inputType}")

            return train, test
        else:
            raise RuntimeError(f"No molecule set attached to model: {repr(self)}. Cannot prepare training data.")

    def getModel(self):
        modelClass = self.trainingStrategy.modelClass
        inputType = self.trainingStrategy.inputType
        if modelClass == "GT" and inputType == "FS":
            model = GraphModel(voc_trg=VocGraph.fromFile(self.vocFile.path))
        elif modelClass == "ST" and inputType == "FS":
            model = GPT2Model(
                VocGPT.fromFile(self.vocFile.path)
            )
        elif modelClass == "SS" and inputType == "MS":
            model = RNN(
                VocSmiles.fromFile(self.vocFile.path),
                is_lstm=True
            )
        else:
            raise NotImplementedError(f"Unknown model configuration: modelClass={modelClass} and inputType={inputType}")

        if self.modelFile:
            model.loadStatesFromFile(self.modelFile.path)

        return model


class DrugExNetValidation(ValidationStrategy):
    validSetSize = models.IntegerField(default=10000)

class DrugExNetTraining(TrainingStrategy):

    class ModelClass(models.TextChoices):
        graphTrans = 'GT', _('Graph Transformer (DrugEx v3)')
        smilesTrans = 'ST', _('SMILES Transformer (DrugEx v3)')
        smilesRNNSingle = 'SS', _('SMILES Single RNN (DrugEx v2)')

    class StructureInputType(models.TextChoices):
        frags = 'FS', _('fragments')
        molecules = 'MS', _('molecules')

    modelClass = models.CharField(max_length=2, choices=ModelClass.choices, default=ModelClass.graphTrans, null=False)
    inputType = models.CharField(max_length=2, choices=StructureInputType.choices, default=StructureInputType.frags, null=False)

class DrugExEnvironment(DataSet):

    class RewardScheme(models.TextChoices):
        paretoCrowding = 'PC', _('Pareto Front with Crowding Distance (PC)')
        paretoSimilarity = 'PS', _('Pareto Front with Similarity (PS)')
        weightedSum = 'WS', _('Weighted Sum (WS)')

    rewardScheme = models.CharField(max_length=2, choices=RewardScheme.choices, default=RewardScheme.paretoCrowding)

    def getInstance(self, use_modifiers=True):
        scorers = []
        thresholds = []
        for scorer in self.scorers.all():
            scorers.append(scorer.getInstance(use_modifiers=use_modifiers))
            thresholds.append(scorer.getThreshold())

        schemes = {
            self.RewardScheme.paretoCrowding: ParetoCrowdingDistance(),
            self.RewardScheme.paretoSimilarity: ParetoSimilarity(),
            self.RewardScheme.weightedSum: WeightedSum()
        }
        reward_scheme = schemes[self.rewardScheme]
        return environment.DrugExEnvironment(scorers, thresholds, reward_scheme)

class DrugExEnvironmentScores(ActivitySet):
    environment = models.ForeignKey(DrugExEnvironment, null=False, on_delete=models.CASCADE, related_name="scores")
    modifiersOn = models.BooleanField(null=False, default=False)

class ScoringMethod(DataSet):

    def getInstance(self, modifier):
        raise NotImplementedError("This method must be overriden in subclasses.")

class ScoreModifier(DataSet):

    def getInstance(self):
        raise NotImplementedError("This method must be overriden in subclasses.")

    @staticmethod
    def test(inputs, **kwargs):
        raise NotImplementedError("This method must be overriden in subclasses.")

class DrugExScorer(DataSet):
    environment = models.ForeignKey(DrugExEnvironment, on_delete=models.CASCADE, null=False, related_name='scorers')
    modifier = models.ForeignKey(ScoreModifier, on_delete=models.CASCADE, null=False)
    method = models.ForeignKey(ScoringMethod, on_delete=models.CASCADE, null=False)
    threshold = models.FloatField(null=False, blank=False)

    def getInstance(self, use_modifiers=True):
        modifier = self.modifier.getInstance() if use_modifiers else None
        return self.method.getInstance(modifier)

    def getThreshold(self):
        return self.threshold

class DrugExAgent(Model):
    environment = models.ForeignKey(DrugExEnvironment, on_delete=models.CASCADE, null=False, related_name='drugexEnviron')
    explorationNet = models.ForeignKey(DrugExNet, on_delete=models.CASCADE, null=False, related_name='drugexExplore') # the prior
    exploitationNet = models.ForeignKey(DrugExNet, on_delete=models.CASCADE, null=False, related_name='drugexExploit') # the agent (the final generator)

    def getGenerator(self):
        return self.generator.all().get() if self.generator.all().exists() else None

    def getDataSetFromMolset(self, molset, data_file, rewrite=True, n_proc=4):
        molset = molset if molset else self.explorationNet.molset
        return self.explorationNet.getDataSetFromMolset(molset, data_file, rewrite, n_proc)


class DrugExAgentValidation(ValidationStrategy):
    validSetSize = models.IntegerField(default=10000)


class DrugExAgentTraining(TrainingStrategy):

    class ExplorerClass(models.TextChoices):
        graph = 'GE', _('Graph Explorer (fragments)')
        smiles = 'SE', _('SMILES Explorer (fragments)')
        smiles_molecules = 'SM', _('SMILES Explorer (molecules)')

    explorer = models.CharField(max_length=2, choices=ExplorerClass.choices, default=ExplorerClass.graph)

    def getExplorerInstance(self, agent, env, mutate, epsilon, beta):
        if self.explorer == self.ExplorerClass.graph:
            return GraphExplorer(agent=agent, env=env, mutate=mutate, epsilon=epsilon, sigma=beta)
        elif self.explorer == self.ExplorerClass.smiles:
            return SmilesExplorer(agent, env, mutate=mutate, epsilon=epsilon, sigma=beta)
        else:
            raise NotImplementedError(f"Unknown explorer class: {self.explorer}")


class DrugEx(Generator):
    agent = models.ForeignKey(Model, on_delete=models.CASCADE, null=False, related_name="generator")
    molset = models.ForeignKey(MolSet, on_delete=models.CASCADE, null=True)
    inputFile = models.FileField(null=True, upload_to='models/', storage=OverwriteStorage())

    def get(self, n_samples):
        import genui.generators.extensions.genuidrugex.genuimodels.builders as builders
        builder_class = getattr(builders, self.agent.builder.name)
        builder = builder_class(Model.objects.get(pk=self.agent.id), noMonitor=True)
        rewrite = False
        if not self.inputFile:
            name = f"DrugEx_{self.pk}_{uuid.uuid4().hex}_input.tsv"
            self.inputFile.save(name, ContentFile(''))
            self.save()
            rewrite = True
        agent = Model.objects.get(pk=self.agent.id)
        data, molset = agent.getDataSetFromMolset(self.molset, self.inputFile, rewrite=rewrite, n_proc=4)
        self.molset = molset
        self.save()
        ret = builder.sample(n_samples, from_inputs=data)[0]
        cleanup()
        return ret

class ModelPerformanceDrugEx(ModelPerfomanceNN):
    isOnValidationSet = models.BooleanField(default=False, blank=False, null=False)
    note = models.CharField(max_length=128, blank=True)

class GenUIModelScorer(ScoringMethod):
    model = models.ForeignKey(QSARModel, on_delete=models.CASCADE, null=False)

    class ModelScorer(Scorer):

        def __init__(self, model, modifier=None):
            super().__init__(modifier=modifier)
            self.model = model
            self.builder = BasicQSARModelBuilder(self.model)

        def getScores(self, mols, frags=None):
            return self.builder.predictMolecules(mols)

        def getKey(self):
            return f"{self.model.name}_{self.model.id}"

    def getInstance(self, modifier):
        return self.ModelScorer(self.model, modifier=modifier)


class PropertyScorer(ScoringMethod):
    _keys = Property(prop='MW').prop_dict.keys()
    _choices = [(key, f'Property: {key}') for key in _keys]

    prop = models.CharField(choices=_choices, max_length=32)

    def getInstance(self, modifier):
        return Property(prop=self.prop, modifier=modifier)

class ClippedScore(ScoreModifier):
    upper = models.FloatField(null=False)
    lower = models.FloatField(null=False, default=0.0)
    high = models.FloatField(null=False, default=1.0)
    low = models.FloatField(null=False, default=0.0)
    smooth = models.BooleanField(null=False, default=False)

    def getInstance(self):
        if self.smooth:
            return modifiers.SmoothClippedScore(self.upper, self.lower, self.high, self.low)
        else:
            return modifiers.ClippedScore(self.upper, self.lower, self.high, self.low)

    @staticmethod
    def test(inputs, **kwargs):
        instance = ClippedScore(**kwargs).getInstance()
        return instance(np.array(inputs))

class SmoothHump(ScoreModifier):
    upper = models.FloatField(null=False, default=1.0)
    lower = models.FloatField(null=False, default=0.0)
    sigma = models.FloatField(null=False, default=0.5)

    def getInstance(self):
        return modifiers.SmoothHump(self.lower, self.upper, self.sigma)

    @staticmethod
    def test(inputs, **kwargs):
        instance = SmoothHump(**kwargs).getInstance()
        return instance(np.array(inputs))