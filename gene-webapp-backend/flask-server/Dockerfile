# Use the Python3.8 image
FROM python:3.8-bullseye

# Set the working directory to /flask
WORKDIR /flask-server

# Copy the app directory contents into the container at /flask
ADD . .

# Install the dependencies
RUN pip install -r requirements.txt

#Expose port
EXPOSE 8000

#ENVIRONMENTAL VARIABLES
ENV FLASK_APP=server.py
ENV FLASK_ENV=production
ENV FLASK_DEBUG=1

CMD ["gunicorn", "server:app", "-w" , "2" , "--threads", "2",  "-b", "0.0.0.0:8000", "--timeout", "600"]


# ENTRYPOINT ["/bin/bash"]

# # run the command to start uWSGI
# CMD ["uwsgi", "--ini", "app.ini"]