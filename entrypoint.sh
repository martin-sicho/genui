#!/usr/bin/env bash

set -e

npm run-script build --prefix /code/frontend
python /code/backend/manage.py migrate --noinput
python /code/backend/manage.py runserver 0.0.0.0:8000