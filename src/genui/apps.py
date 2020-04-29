"""
apps

Created by: Martin Sicho
On: 4/27/20, 6:10 PM
"""
from genui.extensions.utils import discover_extensions

BASE_APPS = [
    'genui.init',
    'genui.projects',
    'genui.compounds',
    'genui.modelling',
    'genui.qsar',
    'genui.generators',
    'genui.maps',
]

def extensions():
    return discover_extensions([
        'genui.extensions',
        'genui.compounds.extensions',
    ])

def all_():
    return BASE_APPS + extensions()

