"""
urls.py in src/genui/search/

"""
from django.urls import path, include
from rest_framework import routers
from . import views


urlpatterns = [
    path('occurrence/inchikey/', views.InchiKeyOccurrenceSearchView.as_view(), name='occurrence_inchikey_search'),
    path('projects/inchikey/', views.InchiKeyProjectsSearchView.as_view(), name='projects_inchikey_search'),

    path('sets/similarity/', views.SimSearchMolsetView.as_view(), name='molsets_sim_search'),
    path('sets/substructure/', views.SubsSearchMolsetView.as_view(), name='molsets_sub_search'),
    path('sets/smarts/', views.SmartsSearchMolsetView.as_view(), name='molsets_smarts_search'),

    path('projects/similarity/', views.SimSearchProjectView.as_view(), name='projects_sim_search'),
    path('projects/substructure/', views.SubsSearchProjectView.as_view(), name='projects_sub_search'),
    path('projects/smarts/', views.SmartsSearchProjectView.as_view(), name='projects_smarts_search'),

    # path('sets/filter/', views.PropertyFilterView.as_view(), name='molsets_property_filters'),
]