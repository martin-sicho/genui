FROM continuumio/miniconda:latest

ARG BASE_DIR="/code/"

ENV PYTHONUNBUFFERED 1
ENV DOCKER_CONTAINER 1

RUN echo "Building the GenUI docker app..."
RUN echo "Current environment:"
RUN printenv

# this is for the celery app -> it needs to wait until the server is up
RUN apt-get install -y --no-install-recommends wait-for-it

COPY ./environment.yml ${BASE_DIR}/environment.yml
RUN conda install python=3.7
RUN conda env update -n base --file ${BASE_DIR}/environment.yml
RUN conda env list
RUN conda list

COPY ./requirements.txt ${BASE_DIR}/requirements.txt
RUN pip install -r ${BASE_DIR}/requirements.txt

COPY ./frontend/package.json ${BASE_DIR}/frontend/package.json
RUN npm --prefix ${BASE_DIR}/frontend install ${BASE_DIR}/frontend

COPY . ${BASE_DIR}/
WORKDIR ${BASE_DIR}/backend/

EXPOSE 8000