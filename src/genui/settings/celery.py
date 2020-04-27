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
celery_app.autodiscover_tasks()