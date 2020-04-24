"""
exceptions

Created by: Martin Sicho
On: 12/27/19, 7:41 PM
"""

class SMILESParsingError(Exception):

    def __init__(self, bad_smiles : str, message : str, *args):
        super().__init__(message, *args)
        self.message = message
        self.bad_smiles = bad_smiles

class StandardizationError(Exception):
    def __init__(self, exp : Exception, message : str, *args):
        super().__init__(message, *args)
        self.message = message
        self.original_exception = exp