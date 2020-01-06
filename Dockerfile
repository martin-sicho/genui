FROM continuumio/miniconda:latest

ENV PYTHONUNBUFFERED 1
ENV DJANGO_ENV dev
ENV DOCKER_CONTAINER 1

RUN conda install python=3.6
RUN conda install -c rdkit rdkit
COPY ./backend/requirements.txt /code/requirements.txt
RUN pip install -r /code/requirements.txt

COPY . /code/
WORKDIR /code/backend/

EXPOSE 8000