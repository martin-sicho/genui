"""
prod

Created by: Martin Sicho
On: 4/27/20, 5:11 PM
"""
import sys

from .base import *

try:
    SECRET_KEY = os.environ['GENUI_BACKEND_SECRET']
except KeyError as exp:
    sys.stderr.write('GENUI_BACKEND_SECRET environment variable must be set in production.\n')
    raise exp

assert type(SECRET_KEY) == str

DEBUG = False

if not GENUI_SETTINGS['HOST_URL']:
    raise KeyError('HOST_URL cannot be empty in production. Use full URI of the host root (i.e. {protocol}://{host}:{port}).')

from .databases.prod import *
