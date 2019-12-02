
# inject the celery app into the genui package
from .celery import celery_app
__all__ = ('celery_app',)
