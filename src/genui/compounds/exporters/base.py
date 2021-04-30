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
        self.listSeparator = ','

    @abstractmethod
    def saveFile(self):
        pass

    def attachActivities(self, rd_mol, mol):
        activity_sets = self.molset.activities
        activities = []
        types = []
        units = []
        sources = []
        sources_id = []
        for activity in mol.activities.all():
            if activity_sets.filter(id=activity.source.id).exists():
                activities.append(str(activity.value))
                types.append(activity.type.value)
                units.append(str(activity.units) if not activity.units else activity.units.value)
                sources.append(activity.source.name)
                sources_id.append(str(activity.source.id))

        rd_mol.SetProp('GENUI_ACTIVITIES', self.listSeparator.join(activities))
        rd_mol.SetProp('GENUI_ACTIVITY_TYPES', self.listSeparator.join(types))
        rd_mol.SetProp('GENUI_ACTIVITY_UNITS', self.listSeparator.join(units))
        rd_mol.SetProp('GENUI_ACTIVITY_SOURCES', self.listSeparator.join(sources))
        rd_mol.SetProp('GENUI_ACTIVITY_SOURCE_IDS', self.listSeparator.join(sources_id))


    def attachMetadata(self, rd_mol, mol):
        rd_mol.SetProp('ID', str(mol.id))
        rd_mol.SetProp('CANONICAL_SMILES', mol.canonicalSMILES)
        rd_mol.SetProp('INCHI', mol.inchi)
        rd_mol.SetProp('INCHI_KEY', mol.inchiKey)

    def getRDMols(self, attachActivities=True, attachMetadata=True):
        ret = []
        for mol, rd_mol in zip(self.molecules, (x.entity.rdMol for x in self.molecules)):
            if attachMetadata:
                self.attachMetadata(rd_mol, mol)
            if attachActivities:
                self.attachActivities(rd_mol, mol)

            ret.append(rd_mol)

        return ret