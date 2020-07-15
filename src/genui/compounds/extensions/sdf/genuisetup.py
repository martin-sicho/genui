"""
genuisetup.py

Created by: Martin Sicho
On: 7/13/20, 1:03 PM
"""
from genui.utils.init import createGroup
from . import models

PARENT = 'genui.compounds'

def setup(*args, **kwargs):
    createGroup(
        "GenUI_Users",
        [
            models.SDFCompounds,
            models.SDFMolecule,
        ],
        force=kwargs['force']
    )

