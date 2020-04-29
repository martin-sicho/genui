from django.apps import AppConfig

class GeneratorsConfig(AppConfig):
    name = 'genui.generators'

    def ready(self, force_inspect=False):
        from . import signals
