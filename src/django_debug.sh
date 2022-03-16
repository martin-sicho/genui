#!/bin/bash

REDIS_NAME=redis_local
POSTGRES_NAME=postgres_local

docker container stop $REDIS_NAME $POSTGRES_NAME
docker container rm $REDIS_NAME $POSTGRES_NAME
docker run --rm -d -p 5432:5432 -v /home/sichom/temp/genui_postgres_data:/var/lib/postgresql/data/ --env POSTGRES_USER=genui --env POSTGRES_PASSWORD=genui --env POSTGRES_DB=genui --name postgres_local informaticsmatters/rdkit-cartridge-debian:Release_2021_03_5 
docker run --rm -d -p 0.0.0.0:6379:6379 --name redis_local redis:alpine 
DJANGO_SETTINGS_MODULE='genui.settings.debug' python manage.py runserver
