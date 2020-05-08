"""
genui

Created by: Martin Sicho
On: 5/7/20, 6:37 PM
"""

import os
import genui.apps

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = 'euws5ei%zq!@0yyo6ta4^e3whylufayu)26th6869x=ljr44=d' if not 'GENUI_BACKEND_SECRET' in os.environ else os.environ['GENUI_BACKEND_SECRET']

# determine if we are running in a docker container
DOCKER = 'DOCKER_CONTAINER' in os.environ and int(os.environ['DOCKER_CONTAINER']) == 1

GENUI_SETTINGS = {
    'HOST' : '',
    'HOST_URL': '',
    'FRONTEND_APP_PATH' : os.environ['GENUI_FRONTEND_APP_PATH'] if 'GENUI_FRONTEND_APP_PATH' in os.environ else None,
    'RF_LOGIN_URL' : 'accounts/rfauth/login/',
    'RF_LOGOUT_URL' : 'accounts/rfauth/logout/',
    'APPS' : genui.apps.all_()
}
# This setting is used to form correct URLs according to the hosting server.
try:
    GENUI_SETTINGS['HOST'] = os.environ['GENUI_BACKEND_HOST']
    GENUI_SETTINGS['HOST_URL'] = f"{os.environ['GENUI_BACKEND_PROTOCOL']}://{GENUI_SETTINGS['HOST']}:{os.environ['GENUI_BACKEND_PORT']}"
except KeyError:
    pass

ALLOWED_HOSTS = []
CSRF_TRUSTED_ORIGINS = []
if GENUI_SETTINGS['HOST']:
    ALLOWED_HOSTS.append(GENUI_SETTINGS['HOST'])
    CSRF_TRUSTED_ORIGINS.append(GENUI_SETTINGS['HOST'])
if DOCKER:
    ALLOWED_HOSTS.append('genui')
    CSRF_TRUSTED_ORIGINS.append('genui')
