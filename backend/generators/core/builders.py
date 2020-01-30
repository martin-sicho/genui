"""
builders

Created by: Martin Sicho
On: 1/26/20, 6:27 PM
"""
import numpy as np
from django.core.files.base import ContentFile
from django.db import transaction

from drugex.api.corpus import CorpusCSV, Corpus
from drugex.api.model.callbacks import PretrainingMonitor
from generators.core.drugex_utils.corpus import CorpusFromDB
from generators.models import DrugeExCorpus
from modelling.core import bases
from generators import models
from . import metrics

class DrugExMonitor(PretrainingMonitor):

    def __init__(self, builder, original_callback=None):
        self.original_call = original_callback
        self.builder = builder
        self.loss_train = None
        self.loss_valid = None
        self.error_rate = None
        self.best_error = None
        self.current_step = None
        self.current_epoch = None
        self.total_steps = None
        self.total_epochs = None
        self.best_state = None
        self.best_yet = False

    def savePerformance(self, metric, value, isValidation, note=""):
        return models.ModelPerformanceDrugEx.objects.create(
            metric=metric,
            value=value,
            isOnValidationSet=isValidation,
            model=self.builder.instance,
            epoch=self.current_epoch,
            step=self.current_step,
            note=note
        )

    @property
    def last_step(self):
        # TODO: determine where the training left off last time

        return 0

    @property
    def last_epoch(self):
        # TODO: determine where the training left off last time

        return 0

    def finalizeStep(self, current_epoch: int, current_batch: int, current_step: int, total_epochs: int,
                     total_batches: int, total_steps: int):

        self.current_step = current_step + self.last_step
        self.current_epoch = current_epoch + self.last_epoch
        self.total_steps = total_steps + self.last_step
        self.total_epochs = total_epochs + self.last_epoch

        info = "Epoch: %d step: %d error_rate: %.3f loss_train: %.3f loss_valid %.3f" % (self.current_epoch, self.current_step, self.error_rate, self.loss_train, self.loss_valid if self.loss_valid is not None else np.inf)
        print(info)

        if self.loss_train is not None:
            self.savePerformance(metrics.DrugExLoss.getDjangoModel(), self.loss_train, False)
        if self.loss_valid is not None:
            self.savePerformance(metrics.DrugExLoss.getDjangoModel(), self.loss_valid, True)
        if self.error_rate:
            if self.best_yet:
                self.savePerformance(metrics.SMILESErrorRate.getDjangoModel(), self.error_rate, False, note=f"Minimum {metrics.SMILESErrorRate.name}, yet.")
            else:
                self.savePerformance(metrics.SMILESErrorRate.getDjangoModel(), self.error_rate, False)

        if self.original_call:
            self.original_call(self)

    def performance(self, loss_train, loss_valid, error_rate, best_error):
        self.loss_train = loss_train
        self.loss_valid = loss_valid
        self.error_rate = error_rate
        self.best_error = best_error

    def smiles(self, smiles, is_valid):
        pass

    def model(self, model):
        pass

    def state(self, current_state, is_best=False):
        if is_best:
            self.best_yet=is_best
            self.best_state = current_state
            self.builder.saveFile()

    def close(self):
        pass

    def getState(self):
        return self.best_state


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
            corpus_init = CorpusCSV.fromFiles(self.initial.corpus.corpusFile.path, self.initial.corpus.vocFile.path)
            voc_all = corpus.voc + corpus_init.voc
            corpus.voc = voc_all
            self.saveCorpusData(corpus)
        return corpus