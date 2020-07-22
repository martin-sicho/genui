"""
genui

Created by: Martin Sicho
On: 5/7/20, 6:37 PM
"""

import os
import genui.apps

# determine if we are running in a docker container
DOCKER = 'DOCKER_CONTAINER' in os.environ and int(os.environ['DOCKER_CONTAINER']) == 1

GENUI_SETTINGS = {
    'HOST' : '',
    'HOST_URL': '',
    'FRONTEND_APP_PATH' : os.environ['GENUI_FRONTEND_APP_PATH'] if 'GENUI_FRONTEND_APP_PATH' in os.environ else None,
    'RF_LOGIN_URL' : 'accounts/rfauth/login/',
    'RF_LOGOUT_URL' : 'accounts/rfauth/logout/',
    'FILES_DIR' : os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../files/')) if 'GENUI_DATA_DIR' not in os.environ else os.path.abspath(os.path.join(os.environ['GENUI_DATA_DIR'], 'files')),
    'APPS' : genui.apps.all_()
}
os.makedirs(GENUI_SETTINGS['FILES_DIR'], exist_ok=True)

# This setting is used to form correct URLs according to the hosting server.
try:
    GENUI_SETTINGS['HOST'] = os.environ['GENUI_BACKEND_HOST']
    GENUI_SETTINGS['HOST_URL'] = f"{os.environ['GENUI_BACKEND_PROTOCOL']}://{GENUI_SETTINGS['HOST']}:{os.environ['GENUI_BACKEND_PORT']}"
except KeyError:
    pass
