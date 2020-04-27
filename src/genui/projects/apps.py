from django.apps import AppConfig

class ProjectsConfig(AppConfig):
    name = 'genui.projects'

    def ready(self, force=False):
        from . import signals
        from . import models
        import django_celery_results.models
        import djcelery_model.models
        from genui.commons.helpers import createGroup

        createGroup(
            "GenUI_Users",
            [
                models.Project
            ],
            force=force
        )

        createGroup(
            "GenUI_Users",
            [
                django_celery_results.models.TaskResult,
                djcelery_model.models.ModelTaskMeta
            ],
            permissions=['view'],
            force=force
        )
