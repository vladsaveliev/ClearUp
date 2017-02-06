# WebServer for fingerprinting
## Install
- `cd /group/ngs/src`
- `git clone /group/ngs/src/Fingerprinting && cd Fingerprinting`
- `virtualenv venv_fp`
- `source venv_fp/bin/activate`
- `pip install requirements.txt`
- `./setup.py develop`

## Set up
- `python manage.py init_db`

## Start server
- `python start.py`

## Add project
- `python manage.py load_project <bcbio_final_path_of_fingerprints_project> <project_name> hg19`

## Load and add project in the AZ US HPC
- `cd /group/ngs/src/Fingerprinting`
- `source /users/klpf990/load_postproc-dev.sh`
- `python manage.py load_project <bcbio_final_path_of_fingerprints_project> <project_name> hg19`
