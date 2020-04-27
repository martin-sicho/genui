from django.apps import AppConfig


class QsarConfig(AppConfig):
    name = 'genui.qsar'

    def ready(self, force_inspect=False):
        from genui.modelling import helpers
        from .core import bases
        helpers.inspectCore('genui.qsar', force=force_inspect, modules=["builders", "descriptors"], additional_bases=[bases.DescriptorCalculator])

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
            force=force_inspect
        )

        createGroup(
            "GenUI_Users",
            [
                models.DescriptorGroup,
            ],
            permissions=['view'],
            force=force_inspect
        )
