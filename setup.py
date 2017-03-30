#!/usr/bin/env python
import os
import pip
from os.path import join, isfile, abspath, dirname, relpath
from setuptools import setup, find_packages

try:
    from ngs_utils import setup_utils
except ImportError:
    pip.main(['install', 'git+git://github.com/vladsaveliev/NGS_Utils.git'])
    from ngs_utils import setup_utils


name = 'Fingerprinting'
package_name = 'fingerprinting'


version = setup_utils.init(name, package_name, __file__)


setup(
    name=name,
    version=version,
    author='Vlad Saveliev, Tsistan Lubinsky, Alla Mikheenko',
    author_email='vladislav.sav@gmail.com',
    description='Mutational fingerprinting tool',
    long_description=(open('README.md').read()),
    keywords='bioinformatics',
    url='https://github.com/AstraZeneca-NGS/Fingerprinting',
    # download_url='https://github.com/AstraZeneca-NGS/Fingerprinting/releases',
    license='GPLv3',
    packages=find_packages(),
    package_data={
        package_name: ['locations/locations.bed'] + setup_utils.find_package_files('prank', package_name),
        '': setup_utils.find_package_files('static', '') + setup_utils.find_package_files('templates', ''),
    },
    include_package_data=True,
    zip_safe=False,
    scripts=['start.py', 'manage.py', 'genotyper.py'],
    install_requires=setup_utils.get_reqs(),
    classifiers=[
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    test_suite='nose.collector',
    tests_require=['nose'],
)

if setup_utils.is_installing():
    print("""
-------------------------
 Installation complete!
-------------------------
Usage:
$ ./genotyper.py bcbio_dir --bed --ref ref --depth 5
$ ./manage.py init_db
$ ./manage.py load_project <bcbio_final_path_of_fingerprints_project>
$ ./start.py
""")
