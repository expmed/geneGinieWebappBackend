# Use the Python3.8 image
FROM python:3.8-bullseye

# Set the working directory to /flask
WORKDIR /flask

# Copy the app directory contents into the container at /flask
ADD . .

# Install the dependencies
RUN pip install -r requirements.txt

# Set the working directory to /app
# WORKDIR /flask/app

RUN chmod +x gunicorn_starter.sh

ENTRYPOINT ["./gunicorn_starter.sh"]

# ENTRYPOINT ["/bin/bash"]

#We can delete this one when used with docker compose
# EXPOSE 8080

# # run the command to start uWSGI
# CMD ["uwsgi", "--ini", "app.ini"]