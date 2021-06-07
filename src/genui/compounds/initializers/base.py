"""
base

Created by: Martin Sicho
On: 18-12-19, 11:32
"""
from abc import ABC, abstractmethod

from django.db import transaction, IntegrityError
from rdkit import Chem

from genui.compounds.initializers.exceptions import SMILESParsingError, StandardizationError, \
    InconsistentIdentifiersException
from genui.compounds.models import MolSet, Molecule, ChemicalEntity
from chembl_structure_pipeline import standardizer as chembl_standardizer

from genui.utils.exceptions import GenUIWarning


class Standardizer(ABC):

    @abstractmethod
    def __call__(self, mol : Chem.Mol):
        pass

class ChEMBLStandardizer(Standardizer):

    def __call__(self, mol):
        if chembl_standardizer.exclude_flag(mol, includeRDKitSanitization=False):
            raise StandardizationError(None, f'ChEMBL standardizer set the exclusion flag for molecule: {Chem.MolToSmiles(mol)}')

        # just for outputs
        try:
            smiles = Chem.MolToSmiles(mol)
        except Exception as exp:
            raise StandardizationError(exp, f'An exception occurred while getting the SMILES for molecule: {mol}')

        try:
            mol = chembl_standardizer.standardize_mol(mol, check_exclusion=False)
        except Exception as exp:
            raise StandardizationError(exp, f'An exception occurred while standardizing molecule: {smiles}')

        try:
            mol, _ = chembl_standardizer.get_parent_mol(mol, check_exclusion=False, verbose=True, neutralize=True)
        except Exception as exp:
            raise StandardizationError(exp, f'An exception occurred while getting the parent molecule of: {smiles}')

        return mol

class MolSetInitializer(ABC):

    class DuplicateChemicalEntityWarning(GenUIWarning):

        def __init__(self, entity, duplicates, compound_set, *args):
            super().__init__(*args)
            self.entity = entity
            self.duplicates = duplicates
            self.compound_set = compound_set

        def getData(self):
            return {
                "duplicates" : [x.id for x in self.duplicates],
                "entity" : {
                    'id' : self.entity.id,
                    'smiles' : self.entity.canonicalSMILES,
                },
                "compound_set" : {
                    'id' : self.compound_set.id,
                    'name' : self.compound_set.name,
                }
            }

    def __init__(self, instance : MolSet, progress_recorder=None, standardizer=None):
        self._instance = instance
        self.standardizer = ChEMBLStandardizer() if not standardizer else standardizer
        self.progress_recorder = progress_recorder
        self.unique_mols = 0
        self.errors = []

    def standardizeFromSMILES(self, smiles):
        rdmol = Chem.MolFromSmiles(smiles, sanitize=True)
        if not rdmol:
            raise SMILESParsingError(smiles, None, f"Failed to create molecule during initialization of molecule set {repr(self._instance)} from SMILES: {smiles}")
        return self.standardizer(rdmol)

    def createMolecule(self, entity, molecule_class, create_kwargs=None):
        if not create_kwargs:
            create_kwargs = dict()

        mol = molecule_class.objects.filter(entity=entity, providers__in=[self.instance])
        if mol.exists():
            try:
                return molecule_class.objects.get(entity=entity, providers__in=[self.instance])
            except molecule_class.MultipleObjectsReturned as exp:
                self.errors.append(self.DuplicateChemicalEntityWarning(entity, mol.all(), self.instance, exp, f"One chemical entity defined for multiple molecules in a compound set."))
                return mol.order_by('id').all()[0] # only return the first instance that was added
        else:
            return molecule_class.objects.create(entity=entity, **create_kwargs)

    def createChemicalEntity(self, smiles):
        rdmol_std = self.standardizeFromSMILES(smiles)
        canon_smiles = Chem.MolToSmiles(rdmol_std, isomericSmiles=True, canonical=True, allHsExplicit=False)
        inchi = Chem.MolToInchi(rdmol_std)
        inchi_key = Chem.InchiToInchiKey(inchi)
        if ChemicalEntity.objects.filter(inchiKey=inchi_key).exists():
            ret = ChemicalEntity.objects.get(
                inchiKey=inchi_key
            )
            print(f'Existing chemical entity found. Returning: {ret.canonicalSMILES} (ID: {ret.id})')
            return ret
        else:
            try:
                print('Creating new chemical entity...')
                if canon_smiles != smiles:
                    print(f'Supplied SMILES string: {smiles} was transformed to a standardized canonical representation: {canon_smiles}')
                return ChemicalEntity.objects.create(
                    canonicalSMILES=canon_smiles,
                    inchi=inchi,
                    inchiKey=inchi_key,
                    rdMol=canon_smiles
                )
            except IntegrityError as exp:
                attempted = {
                    "canonicalSMILES" : canon_smiles,
                    "inchi" : inchi,
                    "inchiKey" : inchi_key,
                }
                existing = {
                    "canonicalSMILES" : '',
                    "inchi" : '',
                    "inchiKey" : ''
                }
                try:
                    mol = ChemicalEntity.objects.get(canonicalSMILES=canon_smiles)
                    existing = {
                        "canonicalSMILES" : mol.canonicalSMILES,
                        "inchi" : mol.inchi,
                        "inchiKey" : mol.inchiKey
                    }
                except Molecule.DoesNotExist:
                    pass
                raise InconsistentIdentifiersException(exp, attempted, existing)

    @transaction.atomic
    def addMoleculeFromSMILES(self, smiles : str, molecule_class=Molecule, create_kwargs=None):
        if not create_kwargs:
            create_kwargs = dict()

        instance = self.getInstance()
        entity = self.createChemicalEntity(smiles)
        molecule = self.createMolecule(entity, molecule_class, create_kwargs)
        molecule.providers.add(instance)
        molecule.save()
        self.unique_mols = instance.molecules.count()

        return molecule

    @abstractmethod
    def populateInstance(self) -> int:
        pass

    @abstractmethod
    def updateInstance(self) -> int:
        pass

    def getInstance(self):
        return self._instance

    @property
    def instance(self):
        return self.getInstance()


