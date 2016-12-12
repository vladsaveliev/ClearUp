# WebServer for fingerprinting

## Install
- `virtualenv venv_fp`
- `source venv_fp/bin/activate`
- `pip install requirements.txt`

## Set up
- `python manage.py init_db`

## Start server
- `python start.py`

## Add project
- `source venv_fp/bin/activate`
- `python manage.py load_project <bcbio_final_path_of_fingerprints_project> <project_name> hg19`

## US HPC
The server is deployed at `/group/ngs/src/Fingerprinting`
