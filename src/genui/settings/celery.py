"""
celery config module

Created by: Martin Sicho
On: 29-11-19, 13:44
"""


import os
from celery import Celery

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'genui.settings.prod')

celery_app = Celery('genui')
celery_app.config_from_object('django.conf:settings', namespace='CELERY')
celery_app.conf.update(
    task_acks_late=True,
    task_track_started=True,
    task_send_sent_event=True,
    worker_prefetch_multiplier=1,
    worker_send_task_events=True,
)
celery_app.autodiscover_tasks()