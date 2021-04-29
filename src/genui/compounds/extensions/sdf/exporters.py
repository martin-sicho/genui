"""
sdf

Created by: Martin Sicho
On: 28.04.21, 13:51
"""

from os.path import join
import time
from django.core.files.base import ContentFile
from rdkit import Chem

from genui.compounds.models import MolSetFile

from genui.compounds.exporters.base import BaseMolSetExporter

class SDFExporter(BaseMolSetExporter):
    name = "SDF File"

    def saveFile(self):
        path = join(f'export_{round(time.time() * 1000)}.sdf')
        file = MolSetFile.objects.create(
            molset=self.molset,
            export=self.instance
        )
        file.file.save(path, ContentFile('placeholder'))

        mols = [x.entity.rdMol for x in self.molecules]
        writer = Chem.SDWriter(open(file.file.path, mode='w'))
        for idx,m in enumerate(mols):
            self.progressRecorder.set_progress(idx, len(mols))
            writer.write(m)
        writer.close()

        file.save()
        return file