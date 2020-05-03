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