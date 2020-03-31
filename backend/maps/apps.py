from django.apps import AppConfig


class MapsConfig(AppConfig):
    name = 'maps'

    def ready(self, force_inspect=False):
        from modelling import helpers
        helpers.inspectCore("maps", force=force_inspect)

        from commons.helpers import createGroup
        from . import models

        createGroup(
            "GenUI_Users",
            [
                models.Map,
                models.MappingStrategy,
                models.Point
            ]
        )
