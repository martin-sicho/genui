"""
genuisetup

Created by: Martin Sicho
On: 7/16/20, 2:48 PM
"""
from genui.utils.init import createGroup
from . import models

PARENT = 'genui.compounds'

def setup(*args, **kwargs):
    createGroup(
        "GenUI_Users",
        [
            models.CSVCompounds,
            models.CSVMolecule,
        ],
        force=kwargs['force']
    )