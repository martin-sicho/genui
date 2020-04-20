#!/usr/bin/env bash

set -e

python /code/backend/manage.py migrate --noinput
npm run-script build --prefix /code/frontend
python /code/backend/manage.py collectstatic --no-input
gunicorn --certfile=/etc/certs/${GENUI_BACKEND_HOST}.crt --keyfile=/etc/certs/${GENUI_BACKEND_HOST}.key genui.wsgi:application --bind 0.0.0.0:${GENUI_BACKEND_PORT}