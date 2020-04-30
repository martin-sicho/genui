"""
genuisetup

Created by: Martin Sicho
On: 4/28/20, 4:38 PM
"""

def setup(*args, **kwargs):
    from genui.utils.init import createGroup
    from . import models

    createGroup(
        "GenUI_Users",
        [
            models.MolSet,
            models.Activity,
            models.ActivitySet,
            models.Molecule,
        ],
        force=kwargs['force']
    )

    createGroup(
        "GenUI_Users",
        [
            models.ActivityTypes,
            models.ActivityUnits,
            models.MoleculePic,
            models.PictureFormat,

        ],
        permissions=['view'],
        force=kwargs['force']
    )

