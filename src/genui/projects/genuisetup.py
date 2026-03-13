
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
