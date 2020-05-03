"""
genuisetup

Created by: Martin Sicho
On: 4/28/20, 4:57 PM
"""

def setup(*args, **kwargs):
    from genui.models import helpers
    from .genuimodels import bases
    helpers.discoverGenuiModels('genui.qsar', force=kwargs['force'], modules=["descriptors"], additional_bases=[bases.DescriptorCalculator])

    from genui.utils.init import createGroup
    from . import models

    createGroup(
        "GenUI_Users",
        [
            models.QSARModel,
            models.ModelActivity,
            models.ModelActivitySet,
            models.QSARTrainingStrategy,
        ],
        force=kwargs['force']
    )

    createGroup(
        "GenUI_Users",
        [
            models.DescriptorGroup,
        ],
        permissions=['view'],
        force=kwargs['force']
    )