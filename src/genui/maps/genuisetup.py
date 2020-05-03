"""
genuisetup

Created by: Martin Sicho
On: 4/28/20, 5:02 PM
"""

def setup(*args, **kwargs):
    from genui.utils.init import createGroup
    from . import models

    createGroup(
        "GenUI_Users",
        [
            models.Map,
            models.MappingStrategy,
            models.Point
        ],
        force=kwargs['force']
    )

