"""
tasks

Created by: Martin Sicho
On: 29-11-19, 13:44
"""

from celery import shared_task

@shared_task(name='addStuff')
def add(a, b):
    return a + b

@shared_task(name='multiplyStuff')
def multiply(a, b):
    return a + b