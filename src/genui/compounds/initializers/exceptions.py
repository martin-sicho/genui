"""
exceptions

Created by: Martin Sicho
On: 12/27/19, 7:41 PM
"""
from genui.utils.exceptions import GenUIException


class SMILESParsingError(GenUIException):

        def __init__(self, bad_smiles : str, original, *args):
            super().__init__(original, *args)
            self.bad_smiles = bad_smiles

        def getData(self):
            return {"bad_smiles" : self.bad_smiles}

class StandardizationError(GenUIException):
    pass

class InconsistentIdentifiersException(GenUIException):

    def __init__(self, original, identifiers, *args):
        super().__init__(original, *args)
        self.identifiers = identifiers

    def getData(self):
        return self.identifiers