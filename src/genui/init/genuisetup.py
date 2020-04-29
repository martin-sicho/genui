"""
genuisetup

Created by: Martin Sicho
On: 4/28/20, 4:33 PM
"""

def setup(*args, **kwargs):
    from genui.commons.helpers import createGroup
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

