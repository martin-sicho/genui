from django.apps import AppConfig


class MapsConfig(AppConfig):
    name = 'maps'

    def ready(self, force_inspect=False):
        from modelling import helpers
        helpers.inspectCore("maps", force=force_inspect)
