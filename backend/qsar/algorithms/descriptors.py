"""
descriptors

Created by: Martin Sicho
On: 16-01-20, 11:08
"""
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
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=bit_len)
                DataStructs.ConvertToNumpyArray(fp, arr)
                fps[i, :] = arr
            except:
                # print(smile) # FIXME: something better should be done in this case
                fps[i, :] = [0] * bit_len
        return pd.DataFrame(fps)

