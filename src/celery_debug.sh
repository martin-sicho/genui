#!/bin/bash

DJANGO_SETTINGS_MODULE='genui.settings.debug' celery -A genui worker --prefetch-multiplier 1 --concurrency 4 -Q celery,gpu --loglevel=info -O fair --hostname celery@%h
