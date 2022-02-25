FROM python:3.8-bullseye

WORKDIR backend

COPY requirements.txt

RUN pip install -r requirements.txt
