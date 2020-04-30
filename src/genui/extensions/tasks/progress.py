"""
progress

Created by: Martin Sicho
On: 4/30/20, 5:29 PM
"""
import json
from _pydecimal import Decimal

import celery_progress.backend
from django.conf import settings
from django.db import IntegrityError
from django_celery_results.models import TaskResult


class ProgressRecorder(celery_progress.backend.ProgressRecorder):

    def set_progress(self, current, total, description="", ctype='application/json', cenc='utf-8'):
        # FIXME: atm some important information is dropped from the original result -> make sure it is transferred to the new one
        try:
            super().set_progress(current, total, description)
        except IntegrityError as exp:
            if hasattr(settings, 'CELERY_TASK_ALWAYS_EAGER') and settings.CELERY_TASK_ALWAYS_EAGER:
                return
            else:
                raise exp
        percent = 0
        if total > 0:
            percent = (Decimal(current) / Decimal(total)) * Decimal(100)
            percent = float(round(percent, 2))
        TaskResult.objects.store_result(
            ctype,
            cenc,
            self.task.request.id,
            json.dumps({
                'current': current,
                'total': total,
                'percent': percent,
                'description': description
            }),
            celery_progress.backend.PROGRESS_STATE,
            task_name=self.task.name,
        )