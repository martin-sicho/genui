"""
signals

Created by: Martin Sicho
On: 1/5/20, 5:40 PM
"""
from django.db.models.signals import pre_delete
from django.dispatch import receiver
from celery import states
from django.conf import settings

from .models import MolSet

@receiver(pre_delete, sender=MolSet, dispatch_uid='on_molset_delete_finish_tasks')
def delete_molset_finish_tasks(sender, instance, using, **kwargs):
    tasks = instance.tasks.all()
    print(tasks)
    for task in tasks:
        task_id = task.task_id
        if task.status not in (states.REVOKED, states.SUCCESS, states.FAILURE):
            settings.CURRENT_CELERY_INSTANCE.control.revoke(task_id=task_id, terminate=True)
    print('About to delete: ', repr(instance))
