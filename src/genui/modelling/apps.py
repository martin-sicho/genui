from django.apps import AppConfig

class ModellingConfig(AppConfig):
    name = 'genui.modelling'

    def ready(self):
        from . import signals
