from django.db import models

# Create your models here.
from django_celery_results.models import TaskResult
from djcelery_model.models import TaskManager
from polymorphic.managers import PolymorphicManager


class TaskShortcutsMixIn:

    def getTasksAsDict(self, started_only=False):
        if started_only:
            tasks = self.tasks.started()
        else:
            tasks = self.tasks.all()

        grouped_tasks = dict()
        for task in tasks:
            task_id = task.task_id
            result = None
            task_name = ''
            try:
                result = TaskResult.objects.get(task_id=task_id)
                task_name = result.task_name
                if not task_name:
                    task_name = ''
            except TaskResult.DoesNotExist:
                task_name = ''
            if task_name not in grouped_tasks:
                grouped_tasks[task_name] = []
            grouped_tasks[task_name].append({
                'task_id' : task_id,
                'status' : task.STATES[task.state][1],
                'result' : result.result if result else None,
                'traceback' : result.traceback if result else None
            })

        return grouped_tasks


class PolymorphicTaskManager(PolymorphicManager, TaskManager):
    pass