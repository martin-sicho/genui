from django.apps import AppConfig

class GeneratorsConfig(AppConfig):
    name = 'generators'

    def ready(self, force_inspect=False):
        from . import signals
        from modelling import helpers
        helpers.inspectCore("generators", force=force_inspect)
