from django.apps import AppConfig

class ProjectsConfig(AppConfig):
    name = 'genui.projects'

    def ready(self, force=False):
        from . import signals
