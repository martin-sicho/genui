"""
builders

Created by: Martin Sicho
On: 1/26/20, 6:27 PM
"""
from django.core.files.base import ContentFile
from django.db import transaction

from drugex.api.corpus import CorpusCSV, Corpus
from generators.core.drugex_utils.corpus import CorpusFromDB
from generators.core.monitors import DrugExNetMonitor
from generators.models import DrugeExCorpus
from modelling.core import bases
from generators import models


class DrugExNetBuilder(bases.ModelBuilder):

    def __init__(self, instance: models.DrugExNet, initial=None, progress=None, onFit=None):
        super().__init__(instance, progress, onFit)
        self.corpus = instance.corpus
        self.initial = initial
        self.onFit = DrugExNetMonitor(self, onFit)

        if not self.corpus:
            self.progressStages.append("Creating Corpus")

    def saveCorpusData(self, corpus : CorpusFromDB):
        prefix = f"{self.algorithmClass.name}{self.instance.id}_project{self.instance.project.id}"
        corpus_name = f"{prefix}_corpus.csv"
        voc_name = f"{prefix}_voc.txt"

        with transaction.atomic():
            self.corpus = DrugeExCorpus.objects.create(network=self.instance)
            self.corpus.vocFile.save(voc_name, ContentFile('placeholder'))
            self.corpus.corpusFile.save(corpus_name, ContentFile('placeholder'))
            corpus.saveVoc(self.corpus.vocFile.path)
            corpus.saveCorpus(self.corpus.corpusFile.path)
            self.corpus.save()

    def getY(self):
        return None

    def getX(self) -> Corpus:
        if not self.corpus:
            self.recordProgress()
            corpus = CorpusFromDB(self.instance.molset)
            corpus.updateData(update_voc=True)
            self.saveCorpusData(corpus)
        else:
            corpus = CorpusCSV.fromFiles(self.corpus.corpusFile.path, self.corpus.vocFile.path)
        if self.initial:
            corpus_init = CorpusCSV.fromFiles(self.initial.corpus.corpusFile.path, self.initial.corpus.vocFile.path)
            voc_all = corpus.voc + corpus_init.voc
            corpus.voc = voc_all
            self.saveCorpusData(corpus)
        return corpus