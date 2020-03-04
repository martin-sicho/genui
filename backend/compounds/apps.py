from django.apps import AppConfig


class CompoundsConfig(AppConfig):
    name = 'compounds'

    def ready(self):
        from . import signals
