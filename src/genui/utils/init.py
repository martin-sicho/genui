"""
init

Created by: Martin Sicho
On: 4/30/20, 7:16 PM
"""
import os
import sys


def checkInitCondition(force):
    if 'GENUI_SKIP_INIT' in os.environ and int(os.environ['GENUI_SKIP_INIT']) == 1:
        return False

    return force or (len(sys.argv) > 1 and sys.argv[1] not in ('makemigrations', 'sqlmigrate', 'migrate'))