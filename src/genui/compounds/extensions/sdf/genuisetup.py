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

