"""
genuisetup

Created by: Martin Sicho
On: 4/28/20, 4:30 PM
"""

def setup(*args, **kwargs):
    from genui.commons.helpers import createGroup
    from . import models

    createGroup(
        "GenUI_Users",
        [
            models.Project
        ],
        force=kwargs['force']
    )