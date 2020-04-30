"""
apps

Created by: Martin Sicho
On: 4/27/20, 6:10 PM
"""
from genui.extensions.utils import discover_extensions

BASE_APPS = [
    'genui.utils',
    'genui.projects',
    'genui.compounds',
    'genui.modelling',
    'genui.qsar',
    'genui.generators',
    'genui.maps',
]

EXTENSION_MODULES = ['genui.extensions'] + [f"{x}.extensions" for x in BASE_APPS]

def extensions():
    return discover_extensions(EXTENSION_MODULES)

def all_():
    return BASE_APPS + extensions()

