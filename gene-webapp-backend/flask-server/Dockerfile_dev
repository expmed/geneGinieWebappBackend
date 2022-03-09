#Base Image
FROM python:3.8-bullseye

WORKDIR /backend

COPY requirements.txt .

RUN pip install -r requirements.txt

COPY ./gene-webapp-backend ./gene-webapp-backend

WORKDIR /backend/gene-webapp-backend

ENV FLASK_APP=server.py

ENV FLASK_ENV=development

ENV export FLASK_DEGUB=1

EXPOSE 5000

CMD ["flask", "run" , "--host=0.0.0.0"]
