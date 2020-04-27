
# inject the celery app into the genui package
from genui.settings.celery import celery_app
__all__ = ('celery_app',)
