#!/usr/bin/env bash

set -e

npm run-script build --prefix /code/frontend
python /code/backend/manage.py migrate --noinput django_rdkit # FIXME: specify a dependency in the migrations of the compounds app rather than this
python /code/backend/manage.py migrate --noinput
python /code/backend/manage.py runserver 0.0.0.0:8000