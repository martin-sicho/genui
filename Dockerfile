FROM continuumio/miniconda:latest

ARG BASE_DIR="/code/"

# setup environment variables
ENV PYTHONUNBUFFERED 1
ENV DOCKER_CONTAINER 1
RUN printenv

# copy over the app files
COPY ./src ${BASE_DIR}/src

# setup dependencies
# this is for the celery worker -> it needs to wait until the server is up
RUN apt-get update && apt-get install -y --no-install-recommends wait-for-it
# initialize the conda environment (we have to use it to fetch rdkit since it is not on pip)
COPY ./environment.yml ${BASE_DIR}/environment.yml
RUN conda install python=3.7 && conda env update -n base --file ${BASE_DIR}/environment.yml && conda env list && conda list
# install the pip packages
COPY ./requirements.txt ${BASE_DIR}/requirements.txt
RUN pip install -r ${BASE_DIR}/requirements.txt
# get the frontend app dependencies
RUN npm --prefix ${BASE_DIR}/src/genui_gui install ${BASE_DIR}/src/genui_gui

# copy the entrypoint script
COPY ./entrypoint.sh ${BASE_DIR}/entrypoint.sh

# set the working directory to where the manage.py lives
WORKDIR ${BASE_DIR}/src