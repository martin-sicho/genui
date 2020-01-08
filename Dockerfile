FROM continuumio/miniconda:latest

ARG BASE_DIR="/code/"

ENV PYTHONUNBUFFERED 1
ENV DOCKER_CONTAINER 1

RUN echo "Building the GenUI docker app..."
RUN echo "Current environment:"
RUN printenv

COPY ./environment.yml ${BASE_DIR}/environment.yml
RUN conda install python=3.7
RUN conda env update -n base --file ${BASE_DIR}/environment.yml
RUN conda env list
RUN conda list

COPY ./requirements.txt ${BASE_DIR}/requirements.txt
RUN pip install -r ${BASE_DIR}/requirements.txt

COPY . ${BASE_DIR}/

WORKDIR ${BASE_DIR}/frontend/
RUN npm install
RUN npm run-script build

WORKDIR ${BASE_DIR}/backend/

EXPOSE 8000