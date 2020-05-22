"""
tasks

Created by: Martin Sicho
On: 5/22/20, 9:09 AM
"""
import time

from celery import shared_task


@shared_task(name='TestTask', bind=True)
def testTask(self):
    for x in range(10):
        print(f"Test Task will sleep: {x}")
        time.sleep(1)
    print("Done.")

@shared_task(name='TestTaskGPU', bind=True, queue='gpu')
def testTaskGPU(self):
    for x in range(10):
        print(f"Test Task on a GPU will sleep: {x}")
        time.sleep(1)
    print("Done.")

