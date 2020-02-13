from django.db import models, transaction

# Create your models here.
from djcelery_model.models import TaskMixin

from commons.models import TaskShortcutsMixIn, PolymorphicTaskManager
from compounds.models import MolSet
from drugex import Voc
from drugex.api.corpus import CorpusCSV, BasicCorpus
from modelling.models import Model, ValidationStrategy, TrainingStrategy, ModelPerfomanceNN, ModelPerformance, ModelFile
from projects.models import DataSet
from qsar.models import QSARModel


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
    def corpus(self):
        corpus_file = self.files.filter(kind=ModelFile.AUXILIARY, note=self.CORPUS_FILE_NOTE)
        voc_file = self.files.filter(kind=ModelFile.AUXILIARY, note=self.VOC_FILE_NOTE)

        if corpus_file:
            corpus_file = corpus_file.get()
        if voc_file:
            voc_file = voc_file.get()

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

class DrugExAgentValidationStrategy(ValidationStrategy):
    pass

class DrugExAgentTrainingStrategy(TrainingStrategy):
    pass

class DrugEx(Generator):
    agent = models.ForeignKey(DrugExAgent, on_delete=models.CASCADE, null=False, related_name="generator")

    def get(self, n_samples):
        from generators.core.builders import DrugExAgentBuilder
        builder = DrugExAgentBuilder(self.agent)
        samples, valids = builder.sample(n_samples)
        return [x for idx, x in enumerate(samples) if bool(valids[idx])]

class ModelPerformanceDrugEx(ModelPerfomanceNN):
    isOnValidationSet = models.BooleanField(default=False, blank=False, null=False)
    note = models.CharField(max_length=128, blank=True)

class ModelPerformanceDrugExAgent(ModelPerformance):
    epoch = models.IntegerField(blank=False, null=False)
    note = models.CharField(max_length=128, blank=True)


