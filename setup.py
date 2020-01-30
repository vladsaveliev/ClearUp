#!/usr/bin/env python
import sys
py_v = sys.version_info[:2]
if not (py_v == (2, 7) or py_v >= (3, 3)):
    sys.exit('Only Python 2.7 or 3.3 and up are supported. Current version: ' + '.'.join(py_v))


package_name = 'clearup'


from ngs_utils import setup_utils
version = setup_utils.init(package_name, package_name, __file__)


from setuptools import setup, find_packages
setup(
    name=package_name,
    version=version,
    author='Vlad Saveliev, Tristan Lubinsky, Alla Miheenko',
    author_email='vladislav.sav@gmail.com',
    description='Mix up check tool',
    long_description=(open('README.md').read()),
    keywords='bioinformatics',
    url='https://github.com/AstraZeneca-NGS/ClearUp',
    license='GPLv3',
    packages=find_packages(),
    package_data={
        package_name:
            ['snps/*.bed.gz*'],
        '':
            setup_utils.find_package_files('static', '') +
            setup_utils.find_package_files('templates', ''),
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

print("""
-------------------------
 Installation complete!
-------------------------
Usage:
$ ./manage.py init_db
$ ./manage.py load_bcbio_project <bcbio_final_path>
$ ./manage.py load_data <dir_with_bams_and_beds> name=<project_name> genome=<genome_build_name>
$ ./start.py
""")
