# Sample mix-up check

## Installation

### Using conda
```bash
conda install -c vladsaveliev -c bioconda -c r -c conda-forge python=3.6 clearup
```

## Set up

- `python manage.py init_db`
- `python manage.py reload_all_pdata`

## Start server

- `python start.py`

## Add project

- `python manage.py load_project <bcbio_final_path_of_fingerprints_project> --name=<project_name>`

## Load and add project in the AZ US HPC

- `cd /gpfs/group/ngs/src/Fingerprinting-1.1`
- `source load.sh`
- `python manage.py load_bcbio_project <bcbio_final_path_of_fingerprints_project> --name=<project_name>`
