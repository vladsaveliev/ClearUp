# WebServer for fingerprinting
## Install
- `pip install --upgrade pip`
- `pip install --upgrade --ignore-installed setuptools`
- `cd /group/ngs/src`
- `git clone /group/ngs/src/Fingerprinting && cd Fingerprinting`
- `virtualenv -p $(which python2) venv_fp`
- `source venv_fp/bin/activate`
- `pip install requirements.txt`
- `./setup.py develop`

## Set up
- `python manage.py init_db`

## Start server
- `python start.py`

## Add project
- `python manage.py load_project <bcbio_final_path_of_fingerprints_project> --name=<project_name>`

## Load and add project in the AZ US HPC
- `cd /group/ngs/src/Fingerprinting`
- `source load.sh`
- `python manage.py load_project <bcbio_final_path_of_fingerprints_project> --name=<project_name>`
