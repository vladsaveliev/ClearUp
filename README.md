# Sample mix-up check

## Install

```bash
source ~/reset_path.sh
module load python
git clone https://github.com/AstraZeneca-NGS/Fingerprinting && cd Fingerprinting
virtualenv -p $(which python2) venv_fp
source venv_fp/bin/activate
pip install --upgrade pip
pip install --upgrade --ignore-installed setuptools
git clone https://github.com/vladsaveliev/NGS_Utils; cd NGS_Utils; pip install -r requirements.txt; ./setup.py develop cd ..
git clone https://github.com/vladsaveliev/NGS_Reporting; cd NGS_Reporting; pip install -r requirements.txt; ./setup.py develop cd ..
git clone https://github.com/AstraZeneca-NGS/reference_data; cd reference_data/bed_annotation; ./setup.py develop; cd ../..
pip install -r requirements.txt
./setup.py develop
module load bcbio-nextgen  #UK and US
echo "" > load.sh
echo "module load python R bedops" >> load.sh  #UK and US
echo "export PATH=`pwd`/venv_fp/bin:`dirname $(dirname $(which bcbio_nextgen.py))`/anaconda/bin:\$PATH" >> load.sh
echo "export PYTHONPATH=" >> load.sh
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
