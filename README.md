# GeneGinie Project

[Project Overview and architecture information](https://drive.google.com/file/d/1LLVk6F5CWMrpBzZxysWPCxlLPgrhoW2i/view?usp=sharing) 

# Basic Setup

## Code Structure

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
    - data
      - citations_db.feather
      - ginie.db
    - requirements.txt
    - server.py
  - nginx
    - Dockerfile
    - nginx.conf
- docker-compose.yml
 ```

## Setup Steps
1. Clone this repo
2. Create a virtual environment (for example)
   2.a python3 -m venv gene-backend
   2.b source gene-backend-env-test/bin/activate
3. Navigate to `./gene-webapp-backend/gene-webapp-backend`
4. Install requirement.txt packages ```pip install -r requirements.txt```
5. navigate to ./flask-server
6. Create a folder named data
7. Copy the database: ginie.db to that folder
8. Copy the dataframe citations_db.feather to that flder
9. `export FLASK_APP=server.py'
10. `flask run -p 8000` We need to use this port because is the one used in the Angular client to create the websocket. If changed it needs to be changed in both sides.  

# Docker Launching for Deploymenmt (TODO)

```
1. DockerCompose Structure
  1.1 Docker Container: flask-server
  1.2 Docker Container: nginx
```
