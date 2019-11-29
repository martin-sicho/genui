"""
tasks

Created by: Martin Sicho
On: 29-11-19, 13:44
"""
from jobs.celery import app

class Point:

    def __init__(self, x, y):
        self.x = x
        self.y = y

@app.task(bind=True, name='addStuff')
def add(self, a, b):
    return Point(a, b)