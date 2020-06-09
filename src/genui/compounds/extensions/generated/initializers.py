"""
generated

Created by: Martin Sicho
On: 07-02-20, 14:55
"""
import traceback

from genui.compounds.initializers.base import MolSetInitializer
from genui.generators.models import Generator
from genui.compounds.extensions.generated.models import GeneratedMolSet


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
        self.nSamples = self.instance.molecules.count()
        self.instance.activities.all().delete()
        self.instance.molecules.clear()
        self.populateInstance()