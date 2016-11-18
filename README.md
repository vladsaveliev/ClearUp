# WebServer for fingerprinting

## Install
- `virtualenv venv_fp`
- `source venv_fp/bin/activate`
- `pip install requirements.txt`

## Set up
- `python manage.py init_db`
- For each project run `python manage.py load_project <bcbio_final_path_of_fingerprints_project>`

## Start server
- `python start.py`