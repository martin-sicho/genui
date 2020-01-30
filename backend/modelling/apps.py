from django.apps import AppConfig

class ModellingConfig(AppConfig):
    name = 'modelling'

    def ready(self, force_inspect=False):
        from . import signals
        from . import helpers
        helpers.inspectCore("modelling", force=force_inspect)
