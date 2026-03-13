from rdkit import Chem
from genui.compounds.extensions.fileimports.parser import FileParser
from . import models

class SDFParser(FileParser):

    def parse(self):
        for mol in Chem.SDMolSupplier(self.path, False, False, False):
            props = {
                'model' : f'{models.__name__}.{models.SDFMolecule.__name__}'
            }
            props.update(mol.GetPropsAsDict())
            props['name'] = mol.GetProp('_Name')
            yield Chem.MolToSmiles(mol), props