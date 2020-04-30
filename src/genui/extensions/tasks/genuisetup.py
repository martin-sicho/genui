"""
genuisetup

Created by: Martin Sicho
On: 4/30/20, 5:51 PM
"""

PARENT = 'genui'

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
