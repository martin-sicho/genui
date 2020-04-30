"""
urls

Created by: Martin Sicho
On: 4/30/20, 5:25 PM
"""
from django.urls import re_path

from . import views

urlpatterns = [
    re_path(r'^api/celery-progress/(?P<task_id>[\w-]+)/$', views.TaskProgressView.as_view())
]
