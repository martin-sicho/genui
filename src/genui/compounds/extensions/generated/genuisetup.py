"""
genuisetup

Created by: Martin Sicho
On: 5/12/20, 9:37 AM
"""

PARENT = 'genui.compounds'

def setup(*args, **kwargs):
    from genui.utils.init import createGroup
    from . import models

    createGroup(
        "GenUI_Users",
        [
            models.GeneratedMolSet,
        ],
        force=kwargs['force']
    )
