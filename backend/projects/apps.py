from django.apps import AppConfig

class ProjectsConfig(AppConfig):
    name = 'projects'

    def ready(self):
        from . import signals
        from . import models
        import django_celery_results.models
        import djcelery_model.models
        from commons.helpers import createGroup

        createGroup(
            "GenUI_Users",
            [
                models.Project
            ]
        )

        createGroup(
            "GenUI_Users",
            [
                django_celery_results.models.TaskResult,
                djcelery_model.models.ModelTaskMeta
            ],
            permissions=['view']
        )
