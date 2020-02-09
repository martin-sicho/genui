"""
generated

Created by: Martin Sicho
On: 07-02-20, 14:55
"""
import time

from compounds.initializers.base import MolSetInitializer
from generators.models import GeneratedMolSet, Generator


class GeneratedSetInitializer(MolSetInitializer):

    def __init__(self, instance: GeneratedMolSet, n_samples, progress=None):
        super().__init__(instance, progress)
        self.nSamples = n_samples

    def populateInstance(self):
        instance = self.getInstance()
        source = Generator.objects.get(pk=instance.source.id)
        smiles = source.get(self.nSamples)
        print(smiles)
        raise NotImplementedError("Implement this first")

    def updateInstance(self):
        if True:
            total = 60
            for i in range(total):
                print(i)
                time.sleep(1)
                if self.progress_recorder:
                    self.progress_recorder.set_progress(i, total)
            return total
        else:
            # FIXME: make this happen
            instance = self.getInstance()
            instance.activities.clear()
            instance.molecules.clear()
            self.populateInstance()