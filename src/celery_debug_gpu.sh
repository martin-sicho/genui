#!/bin/bash

DJANGO_SETTINGS_MODULE='genui.settings.debug' celery -A genui worker --prefetch-multiplier 1 --concurrency 1 -P solo -Q gpu --loglevel=info -O fair --hostname celery@%h
