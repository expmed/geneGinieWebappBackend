#Base Image
FROM python:3.8-bullseye

WORKDIR /backend

COPY requirements.txt .

RUN pip install -r requirements.txt

COPY ./gene-webapp-backend ./gene-webapp-backend

WORKDIR /backend/gene-webapp-backend

CMD ["flask", "run"]
