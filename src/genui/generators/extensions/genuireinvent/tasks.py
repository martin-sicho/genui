from celery import shared_task
import os

from genui.utils.extensions.tasks.progress import ProgressRecorder
from genui.utils.inspection import getObjectAndModuleFromFullName

from django.db import close_old_connections

from . import models
from .torchutils import cleanup


@shared_task(name="BuildReinventModel", bind=True, queue="gpu")
def buildReinventModel(self, model_id, builder_class, model_class):
    # Celery workers should not keep stale DB connections around.
    # Use close_old_connections (not close_all) so tasks can also run eagerly
    # in-process during tests without blowing away the TestCase transaction.
    # close_old_connections()
    try:
        model_cls = getattr(models, model_class)
        instance = model_cls.objects.get(pk=model_id)
        builder_cls = getObjectAndModuleFromFullName(builder_class)[0]
        recorder = ProgressRecorder(self)

        if hasattr(instance, "parent") and instance.parent_id:
            builder = builder_cls(instance, instance.parent, progress=recorder)
        else:
            builder = builder_cls(instance, progress=recorder)

        builder.build()

        return {
            "task_id": getattr(self.request, "id", None),
            "errors": [repr(x) for x in builder.errors],
            "ReinventModelName": instance.name,
            "ReinventModelID": instance.id,
        }
    finally:
        cleanup()
        # close_old_connections()


@shared_task(name="RunReinventStagedLearning", bind=True, queue="gpu")
def runReinventStagedLearning(self, reinvent_id, device="cuda:0"):
    # close_old_connections()
    try:
        instance = models.Reinvent.objects.get(pk=reinvent_id)
        print(f"[TASK] runReinventStagedLearning: Reinvent ID={instance.id}, Project={instance.project_id}")
        print(f"[TASK] Agent ID={instance.agent_id}")

        recorder = ProgressRecorder(self)
        for (cur, desc) in [(0, "Build TOML"), (2, "Finalize")]:
            try:
                recorder.set_progress(cur, 3, description=desc)
            except Exception:
                pass

        print(f"[TASK] Calling instance.run_staged_learning(device={device})")
        toml_path = instance.run_staged_learning(device=device)
        print(f"[TASK] Completed. toml_path={toml_path}")

        rl_log_path = None
        try:
            rl_log_path = instance.agent.get_rl_log_path(generator=instance)
            if rl_log_path and not os.path.exists(rl_log_path):
                rl_log_path = None
        except Exception:
            rl_log_path = None

        try:
            recorder.set_progress(3, 3, description="Done")
        except Exception:
            pass

        return {
            "task_id": getattr(self.request, "id", None),
            "ReinventRunName": instance.name,
            "ReinventRunID": instance.id,
            "toml_path": toml_path,
            "rl_log_path": rl_log_path,
        }
    finally:
        cleanup()
        # close_old_connections()
