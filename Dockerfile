FROM continuumio/miniconda:latest

ARG BASE_DIR="/code/"

# setup environment variables
ENV PYTHONUNBUFFERED 1
ENV DOCKER_CONTAINER 1

# setup dependencies
# wait-for-it: for the worker to wait until the backend is online
# authbind: to allow the app process to bind port 443 if run as nonroot
RUN apt-get update && apt-get install -y --no-install-recommends wait-for-it authbind
# initialize the conda environment (we have to use it to fetch rdkit since it is not on pip)
COPY ./environment.yml ${BASE_DIR}/environment.yml
RUN conda install python=3.7 && conda env update -n base --file ${BASE_DIR}/environment.yml
# install the pip packages
COPY ./requirements.txt ${BASE_DIR}/requirements.txt
RUN pip install -r ${BASE_DIR}/requirements.txt
# get the frontend app dependencies
COPY ./src/genui_gui/package.json ${BASE_DIR}/src/genui_gui/package.json
RUN npm --prefix ${BASE_DIR}/src/genui_gui install ${BASE_DIR}/src/genui_gui

# copy over the sources
COPY ./src ${BASE_DIR}/src

# copy the entrypoint script
COPY ./entrypoint.sh ${BASE_DIR}/entrypoint.sh

# set working directory to where the manage.py lives
WORKDIR ${BASE_DIR}/src