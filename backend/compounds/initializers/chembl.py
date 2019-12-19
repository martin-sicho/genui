"""
chembl

Created by: Martin Sicho
On: 18-12-19, 14:38
"""

from .base import MolSetInitializer

class ChEMBLSetInitializer(MolSetInitializer):

    def populateInstance(self):
        pass

    def getInstance(self):
        return self._instance