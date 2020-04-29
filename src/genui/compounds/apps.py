from django.apps import AppConfig


class CompoundsConfig(AppConfig):
    name = 'genui.compounds'

    def ready(self, force=False):
        from . import signals
