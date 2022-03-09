# Use the Python3.8 image
FROM python:3.8-bullseye

# Set the working directory to /flask
WORKDIR /flask-server

# Copy the app directory contents into the container at /flask
ADD . .

# Install the dependencies
RUN pip install -r requirements.txt

# CMD ["gunicorn server:app -w 2 --threads 2 -b 0.0.0.0:8000"]

# ENTRYPOINT ["/bin/bash"]

#We can delete this one when used with docker compose
EXPOSE 8080