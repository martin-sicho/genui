"""
genuisetup

Created by: Martin Sicho
On: 4/28/20, 4:30 PM
"""

def setup(*args, **kwargs):
    from genui.utils.init import createGroup
    from . import models

    createGroup(
        "GenUI_Users",
        [
            models.Project
        ],
        force=kwargs['force']
    )