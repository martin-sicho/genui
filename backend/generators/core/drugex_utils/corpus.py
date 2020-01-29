"""
corpus

Created by: Martin Sicho
On: 27-01-20, 12:29
"""
from compounds.models import MolSet
from drugex.api.corpus import DataProvidingCorpus
import pandas as pd

class CorpusFromDB(DataProvidingCorpus):

    def __init__(self
                 , molset: MolSet
                 , randomize=True
                 , smiles_field="CANONICAL_SMILES"
                 ):
        super().__init__(smiles_col=smiles_field)
        self.molset = molset
        self.randomize = randomize

    def updateData(self, update_voc=False, sample=None):

        mols = []
        for mol in self.molset.molecules.all():
            mols.append(mol.smiles)

        data = pd.DataFrame(columns=[self.smiles_column])
        data[self.smiles_column] = mols

        # randomize if requested
        if self.randomize:
            data = data.sample(frac=1)

        self.updateDataFrame(data, self.smiles_column, update_voc, sample)

        return self.df, self.voc