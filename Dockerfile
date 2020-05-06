FROM continuumio/miniconda:latest

ARG BASE_DIR="/code/"
ARG GENUI_FRONTEND_APP_BRANCH=master
ARG GENUI_FRONTEND_APP_REPO=git@github.com:martin-sicho/genui_gui.git

# setup environment
ENV PYTHONUNBUFFERED 1
ENV DOCKER_CONTAINER 1
RUN printenv

# this is for the celery app -> it needs to wait until the server is up
RUN apt-get update && apt-get install -y --no-install-recommends wait-for-it
COPY config/.ssh /root/.ssh

COPY ./environment.yml ${BASE_DIR}/environment.yml
RUN conda install python=3.7 && conda env update -n base --file ${BASE_DIR}/environment.yml && conda env list && conda list

COPY ./requirements.txt ${BASE_DIR}/requirements.txt
RUN pip install -r ${BASE_DIR}/requirements.txt

# checkout the frontend app
RUN ssh-keyscan github.com > /root/.ssh/known_hosts && git clone ${GENUI_FRONTEND_APP_REPO} ${BASE_DIR}/src/genui_gui --branch ${GENUI_FRONTEND_APP_BRANCH} && npm --prefix ${BASE_DIR}/src/genui_gui install ${BASE_DIR}/src/genui_gui

COPY ./src/genui ${BASE_DIR}/src/genui
COPY ./src/manage.py ${BASE_DIR}/src/
COPY ./entrypoint.sh ${BASE_DIR}/entrypoint.sh
WORKDIR ${BASE_DIR}/src