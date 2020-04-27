"""
apps

Created by: Martin Sicho
On: 4/27/20, 6:10 PM
"""
from genui.extensions.utils import discover_extensions

BASE_APPS = [
    'genui.projects.apps.ProjectsConfig',
    'genui.compounds.apps.CompoundsConfig',
    'genui.modelling.apps.ModellingConfig',
    'genui.qsar.apps.QsarConfig',
    'genui.generators.apps.GeneratorsConfig',
    'genui.maps.apps.MapsConfig',
]

EXTENSIONS = discover_extensions([
    'genui.extensions',
    'genui.compounds.extensions',
])

ALL = BASE_APPS + EXTENSIONS

