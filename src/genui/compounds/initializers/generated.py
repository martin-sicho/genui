"""
generated

Created by: Martin Sicho
On: 07-02-20, 14:55
"""
import time
import traceback

from compounds.initializers.base import MolSetInitializer
from generators.models import GeneratedMolSet, Generator


class GeneratedSetInitializer(MolSetInitializer):

    def __init__(self, instance: GeneratedMolSet, progress=None, n_samples=100):
        super().__init__(instance, progress)
        self.nSamples = n_samples

    def populateInstance(self):
        instance = self.getInstance()

        if self.progress_recorder:
            self.progress_recorder.set_progress(1, self.nSamples, description="Generating new structures...")
        source = Generator.objects.get(pk=instance.source.id)
        smiles = source.get(self.nSamples)

        for idx, smile in enumerate(smiles):
            try:
                self.addMoleculeFromSMILES(smile)
            except Exception as exp:
                traceback.print_exc()
                self.errors.append(exp)
            if self.progress_recorder:
                self.progress_recorder.set_progress(idx, len(smiles), description=f"Saved {smile}")

        return self.unique_mols

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