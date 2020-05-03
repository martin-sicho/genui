from django.apps import AppConfig

class ModelsConfig(AppConfig):
    name = 'genui.models'

    def ready(self):
        from . import signals
