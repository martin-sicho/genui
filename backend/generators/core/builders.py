"""
builders

Created by: Martin Sicho
On: 1/26/20, 6:27 PM
"""
from django.core.files.base import ContentFile
from django.db import transaction

from drugex.api.corpus import CorpusCSV, Corpus
from drugex.api.model.callbacks import PretrainingMonitor
from generators.core.drugex_utils.corpus import CorpusFromDB
from generators.models import DrugeExCorpus
from modelling.core import bases
from generators import models

class DrugExMonitor(PretrainingMonitor):

    def __init__(self, builder, original_callback=None):
        self.original_call = original_callback
        self.builder = builder
        self.current_model = None

    def finalizeStep(self, current_epoch: int, current_batch: int, current_step: int, total_epochs: int,
                     total_batches: int, total_steps: int):

        if self.original_call:
            self.original_call(self)

    def performance(self, loss_train, loss_valid, error_rate, best_error):
        pass

    def smiles(self, smiles, is_valid):
        pass

    def model(self, model):
        self.current_model = model

    def state(self, current_state, is_best=False):
        if is_best:
            self.builder.saveFile()

    def close(self):
        pass

    def getState(self):
        pass


class DrugExBuilder(bases.ModelBuilder):

    def __init__(self, instance: models.DrugExNet, initial=None, progress=None, onFit=None):
        super().__init__(instance, progress, onFit)
        self.corpus = instance.corpus
        self.initial = initial
        self.onFit = DrugExMonitor(self, onFit)

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
            corpus_init = CorpusCSV.fromFiles(self.corpus.corpusFile.path, self.corpus.vocFile.path)
            voc_all = corpus.voc + corpus_init.voc
            corpus.voc = voc_all
        return corpus