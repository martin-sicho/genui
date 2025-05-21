"""
genuisetup

Created by: Martin Sicho
On: 4/28/20, 9:38 AM
"""

PARENT = 'genui.compounds'

def setup(*args, **kwargs):
    from . import models
    from genui.utils.init import createGroup
    createGroup(
        "GenUI_Users",
        [
            models.PapyrusActivities,
            models.PapyrusActivity,
            models.PapyrusCompounds,
            models.PapyrusMolecule,
        ],
        force=kwargs['force']
    )
    createGroup(
        "GenUI_Users",
        [
            models.PapyrusTarget,
            models.PapyrusAssay,

        ],
        permissions=['view'],
        force=kwargs['force']
    )

