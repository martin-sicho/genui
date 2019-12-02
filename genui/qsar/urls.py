"""
urls

Created by: Martin Sicho
On: 02-12-19, 17:18
"""

from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='index'),
]
