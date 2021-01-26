..  _install-guide-docker:

Deploying a Docker Image
~~~~~~~~~~~~~~~~~~~~~~~~

Likely the easiest way to deploy a production version of the GenUI application
is with the help of Docker images. You can build your own with the GenUI
frontend and backend source code or you can use our Docker images
on Docker Hub. The process of deploying
and building GenUI Docker images is described and implemented in a `separate repository
<https://github.com/martin-sicho/genui-docker>`_ where you will also find a more detailed
overview. However, if you just want to deploy a quick testing
instance of the GenUI app on your local machine, here is a short example:

..  code-block:: bash

    # clone the repository
    git clone https://github.com/martin-sicho/genui-docker
    cd genui-docker

    # get the images
    docker pull sichom/genui-main
    docker pull sichom/genui-gpuworker
    # docker pull sichom/genui-worker # if you do not have an NVIDIA GPU

    # deploy
    GENUI_PROTOCOL=http \
    GENUI_HOST=localhost \
    GENUI_PORT=8000 \
    POSTGRES_PASSWORD=genui \
    REDIS_PASSWORD=redis \
    GENUI_BACKEND_SECRET=`cat django_secret_example` \
    GENUI_CELERY_QUEUES=celery,gpu \
    GENUI_CELERY_CONCURRENCY=`nproc --all` \
    docker-compose -f docker-compose-main.yml -f docker-compose-gpuworker.yml up
    # docker-compose -f docker-compose-main.yml -f docker-compose-worker.yml up # if you do not have an NVIDIA GPU