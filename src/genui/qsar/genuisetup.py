"""
genuisetup

Created by: Martin Sicho
On: 4/28/20, 4:57 PM
"""

def setup(*args, **kwargs):
    from genui.modelling import helpers
    from .core import bases
    helpers.inspectCore('genui.qsar', force=kwargs['force'], modules=["builders", "descriptors"], additional_bases=[bases.DescriptorCalculator])

    from genui.commons.helpers import createGroup
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