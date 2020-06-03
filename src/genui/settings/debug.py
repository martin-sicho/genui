"""
debug.py

Created by: Martin Sicho
On: 4/27/20, 5:10 PM
"""

from .base import *

DEBUG = True

CORS_ORIGIN_WHITELIST = [
"http://localhost:3000", # default frontend host for debugging, change if necessary
]
CORS_ALLOW_CREDENTIALS = True
SESSION_COOKIE_SAMESITE = None

ALLOWED_HOSTS = ['*']

from .databases.debug import *

EMAIL_HOST = 'localhost'
EMAIL_PORT = 1025
EMAIL_HOST_USER = ''
EMAIL_HOST_PASSWORD = ''
EMAIL_USE_TLS = False
DEFAULT_FROM_EMAIL = 'testing@example.com'

MEDIA_ROOT = os.path.join(GENUI_SETTINGS['FILES_DIR'], 'media_debug/')