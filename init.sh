#!/bin/bash

# install docker beforehand
docker run -p 6379:6379 --name genui-redis -d redis
celery worker -A -E genui --loglevel=info &
celery worker -c 2 -A -E genui --loglevel=info --queues gpu &