#Code Structure 

1. DockerCompose
  1.1 Docker Container: flask-server
  1.2 Docker Container: nginx

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
