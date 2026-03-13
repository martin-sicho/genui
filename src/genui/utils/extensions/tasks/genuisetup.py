
PARENT = 'genui.utils'

def setup(*args, **kwargs):
    from genui.utils.init import createGroup
    import django_celery_results.models
    import djcelery_model.models
    createGroup(
        "GenUI_Users",
        [
            django_celery_results.models.TaskResult,
            djcelery_model.models.ModelTaskMeta
        ],
        permissions=['view'],
        force=kwargs['force']
    )
