#!/usr/bin/env bash

set -e

npm run-script build --prefix /code/src/genui_gui/

python manage.py migrate --noinput
python manage.py genuisetup
python manage.py collectstatic --no-input
gunicorn --certfile=/etc/certs/${GENUI_BACKEND_HOST}.crt --keyfile=/etc/certs/${GENUI_BACKEND_HOST}.key genui.wsgi:application --bind 0.0.0.0:${GENUI_BACKEND_PORT}