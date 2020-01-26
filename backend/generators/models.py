from django.db import models

# Create your models here.
from djcelery_model.models import TaskMixin

from commons.models import TaskShortcutsMixIn, PolymorphicTaskManager
from compounds.models import MolSet
from modelling.models import Model, ValidationStrategy
from projects.models import DataSet
from qsar.models import QSARModel


class Generator(TaskShortcutsMixIn, TaskMixin, DataSet):
    objects = PolymorphicTaskManager()

class GeneratedMolSet(MolSet):
    source = models.ForeignKey(Generator, on_delete=models.CASCADE, null=False, related_name="compounds")

class DrugExVocabularyItem(models.Model):
    token = models.CharField(max_length=16, blank=False, null=False)

class DrugExVocabulary(models.Model):
    tokens = models.ManyToManyField(DrugExVocabularyItem)

class DrugeExCorpus(models.Model):
    voc = models.ForeignKey(DrugExVocabulary, on_delete=models.CASCADE, null=False)
    molset = models.ForeignKey(MolSet, on_delete=models.CASCADE, null=False)

class DrugExValidationStrategy(ValidationStrategy):
    validSetSize = models.IntegerField(default=512, null=True)

class DrugExNet(Model):
    corpus = DrugeExCorpus

class DrugExAgent(Model):
    environment = models.ForeignKey(QSARModel, on_delete=models.CASCADE, null=False)
    explorationNet = models.ForeignKey(DrugExNet, on_delete=models.CASCADE, null=True)
    exploitationNet = models.ForeignKey(DrugExNet, on_delete=models.CASCADE, null=True)

class DrugExGenerator(Generator):
    agent = models.ForeignKey(DrugExAgent, on_delete=models.CASCADE, null=False)



