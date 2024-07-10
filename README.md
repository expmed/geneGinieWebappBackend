# Code Structure 

```
1. DockerCompose
  1.1 Docker Container: flask-server
  1.2 Docker Container: nginx
```

```
- gene-webapp-backend
  - env
    - flask-server.env
  - flask-server
    - app
      __init__.py
    - .dockerignore
    - Dockerfile
    - Dockerfile_dev
    - requirements.txt
    - server.py
  - nginx
    - Dockerfile
    - nginx.conf
- docker-compose.yml
- requirements.txt
 ```

# Basic Setup

1. Clone this repo
2. Navigate to `./gene-webapp-backend/gene-webapp-backend`
3. Create a virtual environment
  3.a python3 -m venv gene-backend
  3.b source gene-backend-env-test/bin/activate
4. Install requirement.txt packages ```pip install -r requirements.txt```
5. 

# Docker Launching for Deploymenmt (TODO)
