"""
base

Created by: Martin Sicho
On: 18-12-19, 11:32
"""
from abc import ABC, abstractmethod

from django.db import IntegrityError, transaction
from rdkit import Chem
from rdkit.Chem import AllChem

from genui.compounds.initializers.exceptions import SMILESParsingError, StandardizationError
from genui.compounds.models import MolSet, Molecule
from chembl_structure_pipeline import standardizer as chembl_standardizer


class Standardizer(ABC):

    @abstractmethod
    def __call__(self, mol : Chem.Mol):
        pass

class ChEMBLStandardizer(Standardizer):

    def __call__(self, mol):
        if chembl_standardizer.exclude_flag(mol, includeRDKitSanitization=False):
            raise StandardizationError(None, f'ChEMBL standardizer set the exclusion flag for molecule: {Chem.MolToSmiles(mol)}')

        try:
            smiles = Chem.MolToSmiles(mol)
        except Exception as exp:
            raise StandardizationError(exp, f'An exception occurred while getting the SMILES for molecule: {mol}')

        try:
            mol = chembl_standardizer.standardize_mol(mol, check_exclusion=False)
        except Exception as exp:
            raise StandardizationError(exp, f'An exception occurred while standardizing: {smiles}')

        try:
            mol, _ = chembl_standardizer.get_parent_mol(mol, check_exclusion=False, verbose=True, neutralize=True)
        except Exception as exp:
            raise StandardizationError(exp, f'An exception occurred while getting the parent molecule of: {smiles}')

        return mol

class MolSetInitializer(ABC):

    def __init__(self, instance : MolSet, progress_recorder=None, standardizer=None):
        self._instance = instance
        self.standardizer = ChEMBLStandardizer() if not standardizer else standardizer
        self.progress_recorder = progress_recorder
        self.unique_mols = 0
        self.errors = []

    def addMoleculeFromSMILES(self, smiles : str, molecule_class=Molecule, constructor_kwargs=None):
        # TODO: check if molecule_class is a subclass of Molecule
        if not constructor_kwargs:
            constructor_kwargs = dict()

        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if not mol:
            raise SMILESParsingError(smiles, f"Failed to create molecule during initialization of molecule set {repr(self._instance)} from SMILES: {smiles}")
        smol = self.standardizer(mol)
        canon_smiles = Chem.MolToSmiles(smol, isomericSmiles=True, canonical=True, allHsExplicit=True)
        inchi_key = Chem.MolToInchiKey(smol)
        params = {
            "canonicalSMILES" : canon_smiles
            , "inchiKey" : inchi_key
        }
        params.update(constructor_kwargs)

        with transaction.atomic():
            instance = self.getInstance()
            if molecule_class.objects.filter(**params).count():
                ret = molecule_class.objects.get(**params)
                # TODO: add a callback or something to check if there are no issues (like two CHMEBL molecules that are the same but have different ID)
            else:
                try:
                    params["molObject"] = smol
                    ret = molecule_class.objects.create(**params)
                    ret.providers.add(self._instance)
                    ret.morganFP2 = AllChem.GetMorganFingerprintAsBitVect(ret.molObject, radius=2, nBits=512)
                    ret.save()
                except IntegrityError as exp:
                    # TODO: analyze the error and provide more details to the caller
                    raise exp

            ret.providers.add(instance)
            ret.save()
        self.unique_mols = instance.molecules.count()
        return ret

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


