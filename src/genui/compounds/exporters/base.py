"""
base

Created by: Martin Sicho
On: 28.04.21, 13:45
"""
from abc import ABC, abstractmethod


class BaseMolSetExporter(ABC):
    name = None

    def __init__(self, export_instance, progress_recorder=None):
        self.instance = export_instance
        self.molset = export_instance.molset
        self.molecules = self.molset.molecules.all()
        self.progressRecorder = progress_recorder
        self.errors = []

    @abstractmethod
    def saveFile(self):
        pass