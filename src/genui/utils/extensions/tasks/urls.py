from django.urls import re_path

from . import views

urlpatterns = [
    re_path(r'^api/tasks/progress/(?P<task_id>[\w-]+)/$', views.TaskProgressView.as_view())
]
