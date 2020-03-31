"""
builders

Created by: Martin Sicho
On: 1/26/20, 6:27 PM
"""
from django.core.files.base import ContentFile
from django.db import transaction
from pandas import Series

from drugex.api.corpus import CorpusCSV, Corpus, BasicCorpus
from generators.core.drugex_utils.corpus import CorpusFromDB
from generators.core.monitors import DrugExNetMonitor, DrugExAgentMonitor
from modelling.core import bases
from generators import models
from modelling.models import ModelFile


class DrugExNetBuilder(bases.ProgressMixIn, bases.ModelBuilder):

    def __init__(self, instance: models.DrugExNet, initial: models.DrugExNet=None, progress=None, onFit=None):
        super().__init__(instance, progress, onFit)
        self.corpus = instance.corpus
        self.initial = initial
        self.onFit = DrugExNetMonitor(self, onFit)

        if not self.corpus:
            self.progressStages.append("Creating Corpus")

    def createCorpus(self):
        if self.instance.molset:
            corpus = CorpusFromDB(self.instance.molset)
            corpus.updateData(update_voc=True)
            with transaction.atomic():
                if self.instance.corpus:
                    self.instance.corpus = None
                corpus_file = ModelFile.create(
                    self.instance,
                    "corpus.csv",
                    ContentFile('placeholder'),
                    note=models.DrugExNet.CORPUS_FILE_NOTE
                )
                voc_file = ModelFile.create(
                    self.instance,
                    "voc.txt",
                    ContentFile('placeholder'),
                    note=models.DrugExNet.VOC_FILE_NOTE
                )
                corpus.saveVoc(voc_file.path)
                corpus.saveCorpus(corpus_file.path)
                self.corpus = self.instance.corpus
        elif self.instance.corpus:
            print("WARNING:  No molset available to create corpus. Falling back to the original...")
            self.corpus = self.instance.corpus
        else:
            Exception("Unable to create corpus. No molecule set is specified and no corpus found on model instance.")

    def getY(self):
        return None

    def getX(self) -> Corpus:
        if not self.corpus:
            self.recordProgress()
            self.createCorpus()

        if self.initial:
            corpus_init = self.initial.corpus
            # voc_all = self.corpus.voc + corpus_init.voc
            # self.corpus.voc = voc_all
            # FIXME: add an error message if there are extra tokens in this vocabulary when compared to parent
            self.corpus.voc = corpus_init.voc
            self.corpus.saveVoc(self.instance.vocFile.path)

        return self.corpus

    def build(self) -> models.Model:
        if self.instance.molset and self.validation:
            return super().build()
        else:
            raise NotImplementedError("Building DrugEx network without molecule set and validation strategy is not allowed.")

    def sample(self, n_samples):
        return self.model.sample(n_samples)

class DrugExAgentBuilder(bases.ProgressMixIn, bases.ModelBuilder):

    def __init__(
            self,
            instance: models.DrugExAgent,
            progress=None,
            onFit=None
    ):
        super().__init__(instance, progress, onFit)
        self.onFit = DrugExAgentMonitor(self, self.onFit)
        self.exploitNet = self.instance.exploitationNet
        self.exploreNet = self.instance.explorationNet
        self.environ = self.instance.environment
        self.corpus = BasicCorpus(vocabulary=self.exploitNet.corpus.voc)

    def getY(self) -> Series:
        pass

    def getX(self) -> Corpus:
        return self.corpus

    def sample(self, n_samples):
        return self.model.sample(n_samples)