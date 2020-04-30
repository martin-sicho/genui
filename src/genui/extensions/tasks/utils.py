"""
tasks

Created by: Martin Sicho
On: 1/4/20, 7:25 PM
"""


def runTask(task, instance=None, eager=False, args=tuple(), kwargs=dict()):
    if eager:
        result = task(*args, **kwargs)
    elif instance:
        result = instance.apply_async(task, args=args, kwargs=kwargs)
    else:
        result = task.apply_async(args=args, kwargs=kwargs)

    task_id = result.id if hasattr(result, 'id') else None
    return result, task_id
