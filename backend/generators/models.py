from django.db import models

# Create your models here.
from djcelery_model.models import TaskMixin

from commons.models import TaskShortcutsMixIn, PolymorphicTaskManager
from compounds.models import MolSet
from modelling.models import Model, ValidationStrategy, TrainingStrategy
from projects.models import DataSet
from qsar.models import QSARModel


class Generator(TaskShortcutsMixIn, TaskMixin, DataSet):
    objects = PolymorphicTaskManager()

class GeneratedMolSet(MolSet):
    source = models.ForeignKey(Generator, on_delete=models.CASCADE, null=False, related_name="compounds")

# class DrugExVocabularyItem(models.Model):
#     token = models.CharField(max_length=16, blank=False, null=False)
#
# class DrugExVocabulary(models.Model):
#     tokens = models.ManyToManyField(DrugExVocabularyItem)

class DrugeExCorpus(models.Model):
    # voc = models.ForeignKey(DrugExVocabulary, on_delete=models.CASCADE, null=False)
    corpusFile = models.FileField(null=True, blank=True, upload_to='drugex/corpus/')
    vocFile = models.FileField(null=True, blank=True, upload_to='drugex/corpus/')

class DrugExNet(Model):
    molset = models.ForeignKey(MolSet, on_delete=models.CASCADE, null=False)
    corpus = models.ForeignKey(DrugeExCorpus, on_delete=models.CASCADE, null=True)

class DrugExValidationStrategy(ValidationStrategy):
    validSetSize = models.IntegerField(default=512, null=True)

class DrugExTrainingStrategy(TrainingStrategy):
    pass

class DrugExGenerator(Generator):
    environment = models.ForeignKey(QSARModel, on_delete=models.CASCADE, null=False, related_name='drugexEnviron')
    explorationNet = models.ForeignKey(DrugExNet, on_delete=models.CASCADE, null=False, related_name='drugexExplore')
    exploitationNet = models.ForeignKey(DrugExNet, on_delete=models.CASCADE, null=False, related_name='drugexExploit')



