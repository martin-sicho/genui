"""
base

Created by: Martin Sicho
On: 18-12-19, 11:32
"""
from abc import ABC, abstractmethod

from compounds.models import MolSet


class MolSetInitializer(ABC):

    def __init__(self, instance : MolSet):
        self._instance = instance

    @abstractmethod
    def populateInstance(self):
        pass

    @abstractmethod
    def getInstance(self):
        pass


