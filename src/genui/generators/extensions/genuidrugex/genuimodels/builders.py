"""
builders

Created by: Martin Sicho
On: 1/26/20, 6:27 PM
"""
from django.core.files.base import ContentFile
from django.db import transaction
from pandas import Series
import torch

from drugex.api.corpus import Corpus, BasicCorpus
from genui.utils import gpu
from genui.generators.extensions.genuidrugex.genuimodels.corpus import CorpusFromDB
from .monitors import DrugExNetMonitor, DrugExAgentMonitor
from genui.models.genuimodels import bases
from genui.models.models import ModelFile, Model
from ..models import DrugExNet, DrugExAgent
from ..torchutils import cleanup

class DrugExBuilderMixIn:

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.device = None
        self.allocateDevice()

    def __del__(self):
        if self.device:
            self.device = None
            cleanup()
            self.releaseDevice()

    def releaseDevice(self):
        if self.device and self.device != 'cpu':
            print(f'Releasing device: {self.device}')
            gpu.release(self.device)
            self.device = None
        else:
            self.device = None

    def allocateDevice(self):
        if not self.device:
            self.device = gpu.allocate() # TODO: wait for some time and try again if we get an allocation exception
            if not self.device:
                print('Failed to allocate GPU device. Using CPU...')
                self.device = 'cpu'
                torch.device(self.device)
            else:
                torch.device('cuda', int(self.device['index']))

class DrugExNetBuilder(bases.ProgressMixIn, DrugExBuilderMixIn, bases.ModelBuilder):

    def __init__(self, instance: DrugExNet, initial: DrugExNet =None, progress=None, onFit=None):
        super().__init__(instance, progress, onFit)
        self.corpus = instance.corpus
        self.initial = initial
        self.onFit = DrugExNetMonitor(self, onFit)

        if not self.corpus:
            self.progressStages.append("Creating Corpus")

    @property
    def corePackage(self):
        from .. import genuimodels
        return genuimodels

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
                    note=DrugExNet.CORPUS_FILE_NOTE
                )
                voc_file = ModelFile.create(
                    self.instance,
                    "voc.txt",
                    ContentFile('placeholder'),
                    note=DrugExNet.VOC_FILE_NOTE
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

    def build(self) -> Model:
        if self.instance.molset and self.validation:
            return super().build()
        else:
            raise NotImplementedError("Building DrugEx network without molecule set and validation strategy is not allowed.")

    def sample(self, n_samples):
        return self.model.sample(n_samples)

class DrugExAgentBuilder(bases.ProgressMixIn, DrugExBuilderMixIn, bases.ModelBuilder):

    def __init__(
            self,
            instance: DrugExAgent,
            progress=None,
            onFit=None
    ):
        super().__init__(instance, progress, onFit)
        self.onFit = DrugExAgentMonitor(self, self.onFit)
        self.exploitNet = self.instance.exploitationNet
        self.exploreNet = self.instance.explorationNet
        self.environ = self.instance.environment
        self.corpus = BasicCorpus(vocabulary=self.exploitNet.corpus.voc)

    @property
    def corePackage(self):
        from .. import genuimodels
        return genuimodels

    def getY(self) -> Series:
        pass

    def getX(self) -> Corpus:
        return self.corpus

    def sample(self, n_samples):
        return self.model.sample(n_samples)