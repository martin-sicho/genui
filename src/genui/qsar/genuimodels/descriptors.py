"""
descriptors

Created by: Martin Sicho
On: 16-01-20, 11:08
"""
import traceback

from rdkit.Chem.Scaffolds import MurckoScaffold

from . import bases
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import pandas as pd


class MorganFPCalculator(bases.DescriptorCalculator):
    group_name = "MORGANFP"

    def __call__(self, smiles, radius=3, bit_len=4096, scaffold=0):
        fps = np.zeros((len(smiles), bit_len))
        for i, smile in enumerate(smiles):
            mol = Chem.MolFromSmiles(smile)
            arr = np.zeros((1,))
            try:
                if scaffold == 1:
                    mol = MurckoScaffold.GetScaffoldForMol(mol)
                elif scaffold == 2:
                    mol = MurckoScaffold.MakeScaffoldGeneric(mol)
                if not mol:
                    raise Exception(f'Failed to calculate Morgan fingerprint (creating RDKit instance from smiles failed: {smile})')
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=bit_len)
                DataStructs.ConvertToNumpyArray(fp, arr)
                fps[i, :] = arr
            except Exception as exp:
                # TODO: use a more specific exception related to descriptor errors
                # traceback.print_exc()
                self.builder.errors.append(exp)
                fps[i, :] = [0] * bit_len
        return pd.DataFrame(fps)

