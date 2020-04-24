from django.apps import AppConfig


class QsarConfig(AppConfig):
    name = 'qsar'

    def ready(self, force_inspect=False):
        from modelling import helpers
        from .core import bases
        helpers.inspectCore("qsar", force=force_inspect, modules=["builders", "descriptors"], additional_bases=[bases.DescriptorCalculator])

        from commons.helpers import createGroup
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
