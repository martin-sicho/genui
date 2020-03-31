#!/usr/bin/env bash

set -e

python /code/backend/manage.py migrate --noinput
npm run-script build --prefix /code/frontend
python /code/backend/manage.py runserver 0.0.0.0:8000