#!/bin/bash

# install docker beforhand
docker run -p 6379:6379 --name genui-redis -d redis
celery -A jobs.celery worker --loglevel=info
