VENV_NAME=venv_fp

source ~/reset_path.sh  #UK and US
module load python  #UK and US
virtualenv ${VENV_NAME}
source ${VENV_NAME}/bin/activate
pip install --upgrade pip
pip install --upgrade --ignore-installed setuptools
git clone --recursive https://github.com/vladsaveliev/NGS_Utils; cd NGS_Utils; pip install -r requirements.txt; ./setup.py develop; cd ..
git clone --recursive https://github.com/AstraZeneca-NGS/reference_data; cd reference_data/bed_annotation; ./setup.py develop; cd ../..
git clone --recursive https://github.com/AstraZeneca-NGS/NGS_Reporting ; cd NGS_Reporting;  pip install -r requirements.txt; ./setup.py develop; cd ..
pip install -r requirements.txt; ./setup.py develop
module load bcbio-nextgen  #UK and US
echo "module load python" > load.sh  #UK and US
echo "export PATH=`pwd`/${VENV_NAME}/bin:`dirname $(which bcbio_nextgen.py)`:\$PATH" >> load.sh
echo "export PYTHONPATH=" >> load.sh
