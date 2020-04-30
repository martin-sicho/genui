from django.db import models, transaction

# Create your models here.
from djcelery_model.models import TaskMixin

from genui.extensions.tasks.models import TaskShortcutsMixIn, PolymorphicTaskManager
from genui.compounds.models import MolSet
from drugex import Voc
from drugex.api.corpus import CorpusCSV, BasicCorpus
from genui.modelling.models import Model, ValidationStrategy, TrainingStrategy, ModelPerfomanceNN, ModelPerformance, ModelFile
from genui.projects.models import DataSet
from genui.qsar.models import QSARModel


class Generator(TaskShortcutsMixIn, TaskMixin, DataSet):
    objects = PolymorphicTaskManager()

    # TODO: it would be useful to have this as an abstract method if possible
    def get(self, n_samples) -> [str]:
        """
        All generators should override this. This should
        return a list of SMILES strings of size N.

        :return: list of SMILES strings
        """

        raise NotImplemented("You have to override this method in subclasses.")

class GeneratedMolSet(MolSet):
    source = models.ForeignKey(Generator, on_delete=models.CASCADE, null=False, related_name="compounds")

class DrugExNet(Model):
    CORPUS_FILE_NOTE = "drugex_corpus"
    VOC_FILE_NOTE = "drugex_voc"

    molset = models.ForeignKey(MolSet, on_delete=models.CASCADE, null=True)
    parent = models.ForeignKey("self", on_delete=models.CASCADE, null=True)

    @property
    def corpusFile(self):
        ret = self.files.filter(kind=ModelFile.AUXILIARY, note=self.CORPUS_FILE_NOTE)
        return ret.get() if ret else None

    @property
    def vocFile(self):
        ret = self.files.filter(kind=ModelFile.AUXILIARY, note=self.VOC_FILE_NOTE)
        return ret.get() if ret else None

    @property
    def corpus(self):
        corpus_file = self.corpusFile
        voc_file = self.vocFile

        corpus = None
        if voc_file and corpus_file:
            corpus = CorpusCSV.fromFiles(corpus_file.path, voc_file.path)
        elif voc_file:
            corpus = BasicCorpus(vocabulary=Voc(voc_file.path))

        return corpus

# class DrugExVocabularyItem(models.Model):
#     token = models.CharField(max_length=16, blank=False, null=False)
#
# class DrugExVocabulary(models.Model):
#     tokens = models.ManyToManyField(DrugExVocabularyItem)

class DrugExValidationStrategy(ValidationStrategy):
    validSetSize = models.IntegerField(default=512, null=True)

class DrugExNetTrainingStrategy(TrainingStrategy):
    pass

class DrugExAgent(Model):
    environment = models.ForeignKey(QSARModel, on_delete=models.CASCADE, null=False, related_name='drugexEnviron')
    explorationNet = models.ForeignKey(DrugExNet, on_delete=models.CASCADE, null=False, related_name='drugexExplore')
    exploitationNet = models.ForeignKey(DrugExNet, on_delete=models.CASCADE, null=False, related_name='drugexExploit')

    def getGenerator(self):
        return self.generator.all().get() if self.generator.all().exists() else None

class DrugExAgentValidationStrategy(ValidationStrategy):
    pass

class DrugExAgentTrainingStrategy(TrainingStrategy):
    pass

class DrugEx(Generator):
    agent = models.ForeignKey(Model, on_delete=models.CASCADE, null=False, related_name="generator")

    def get(self, n_samples):
        import genui.generators.core.builders as builders
        builder_class = getattr(builders, self.agent.builder.name)
        builder = builder_class(Model.objects.get(pk=self.agent.id))
        samples, valids = builder.sample(n_samples)
        return [x for idx, x in enumerate(samples) if bool(valids[idx])]

class ModelPerformanceDrugEx(ModelPerfomanceNN):
    isOnValidationSet = models.BooleanField(default=False, blank=False, null=False)
    note = models.CharField(max_length=128, blank=True)

class ModelPerformanceDrugExAgent(ModelPerformance):
    epoch = models.IntegerField(blank=False, null=False)
    note = models.CharField(max_length=128, blank=True)


