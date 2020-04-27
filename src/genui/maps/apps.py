from django.apps import AppConfig


class MapsConfig(AppConfig):
    name = 'genui.maps'

    def ready(self, force_inspect=False):
        from genui.modelling import helpers
        helpers.inspectCore('genui.maps', force=force_inspect)

        from genui.commons.helpers import createGroup
        from . import models

        createGroup(
            "GenUI_Users",
            [
                models.Map,
                models.MappingStrategy,
                models.Point
            ],
            force=force_inspect
        )
