"""
worker

Created by: Martin Sicho
On: 29-11-19, 13:44
"""


from celery import Celery

app = Celery(
  broker='redis://localhost',
  backend='redis://localhost',
  include=['jobs.tasks']
)


