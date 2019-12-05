"""
urls.py

Created by: Martin Sicho
On: 04-12-19, 15:06
"""

from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='projects-index'),
]