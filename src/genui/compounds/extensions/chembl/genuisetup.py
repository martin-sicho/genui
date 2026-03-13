PARENT = 'genui.compounds'

def setup(*args, **kwargs):
    from . import models
    from genui.utils.init import createGroup
    createGroup(
        "GenUI_Users",
        [
            models.ChEMBLActivities,
            models.ChEMBLActivity,
            models.ChEMBLCompounds,
            models.ChEMBLMolecule,
        ],
        force=kwargs['force']
    )
    createGroup(
        "GenUI_Users",
        [
            models.ChEMBLTarget,
            models.ChEMBLAssay,

        ],
        permissions=['view'],
        force=kwargs['force']
    )

