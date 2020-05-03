"""
apps

Created by: Martin Sicho
On: 4/27/20, 6:10 PM
"""
from genui.utils.inspection import discover_extensions

BASE_APPS = [
    'genui.utils',
    'genui.accounts',
]

API_APPS = [
    'genui.projects',
    'genui.compounds',
    'genui.modelling',
    'genui.qsar',
    'genui.generators',
    'genui.maps',
]

def extensions():
    return discover_extensions(['genui.extensions'] + [f"{x}.extensions" for x in BASE_APPS + API_APPS])

def all_():
    return BASE_APPS + API_APPS + extensions()

