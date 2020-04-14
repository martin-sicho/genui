"""
models

Created by: Martin Sicho
On: 1/12/20, 3:16 PM
"""
import os

from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.db import models
from django_celery_results.models import TaskResult
from djcelery_model.models import TaskManager
from polymorphic.managers import PolymorphicManager

class OverwriteStorage(FileSystemStorage):

    def get_available_name(self, name, **kwargs):
        """Returns a filename that's free on the target storage system, and
        available for new content to be written to.

        Found at http://djangosnippets.org/snippets/976/

        This file storage solves overwrite on upload problem. Another
        proposed solution was to override the save method on the model
        like so (from https://code.djangoproject.com/ticket/11663):

        def save(self, *args, **kwargs):
            try:
                this = MyModelName.objects.get(id=self.id)
                if this.MyImageFieldName != self.MyImageFieldName:
                    this.MyImageFieldName.delete()
            except: pass
            super(MyModelName, self).save(*args, **kwargs)
        """
        # If the filename already exists, remove it as if it was a true file system
        if self.exists(name):
            os.remove(os.path.join(settings.MEDIA_ROOT, name))
        return name

class TaskShortcutsMixIn:

    def getTasksAsDict(self, started_only=False):
        if started_only:
            tasks = self.tasks.started()
        else:
            tasks = self.tasks.all()

        grouped_tasks = dict()
        for task in tasks:
            task_id = task.task_id
            try:
                result = TaskResult.objects.get(task_id=task_id)
                task_name = result.task_name
            except TaskResult.DoesNotExist:
                print(f"Task {task_id} not found in the database. Skipping...")
                continue
            if not task_name:
                task_name = 'UnknownTask'
            if task_name not in grouped_tasks:
                grouped_tasks[task_name] = []
            grouped_tasks[task_name].append(result)

        data = dict()
        for key in grouped_tasks:
            data[key] = [{"task_id" : x.task_id, "status" : x.status, "result": x.result, "traceback": x.traceback} for x in grouped_tasks[key]]
        return data


class PolymorphicTaskManager(PolymorphicManager, TaskManager):
    pass

def NON_POLYMORPHIC_CASCADE(collector, field, sub_objs, using):
    """
    This is a special cascade implementation to fix some delete errors
    when cascading polymorphic models.

    See: https://github.com/django-polymorphic/django-polymorphic/issues/229#issuecomment-398434412

    :param collector:
    :param field:
    :param sub_objs:
    :param using:
    :return:
    """

    return models.CASCADE(collector, field, sub_objs.non_polymorphic(), using)