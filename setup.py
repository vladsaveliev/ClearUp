#!/usr/bin/env python
import sys
py_v = sys.version_info[:2]
if not (py_v == (2, 7) or py_v >= (3, 3)):
    sys.exit('Only Python 2.7 or 3.3 and up are supported. Current version: ' + '.'.join(py_v))
    

import clearup
package_name = clearup.__name__


from ngs_utils import setup_utils
version = setup_utils.init(package_name, package_name, __file__)


from setuptools import setup, find_packages
setup(
    name=package_name,
    version=version,
    author='Vlad Saveliev, Tsistan Lubinsky, Alla Mikheenko',
    author_email='vladislav.sav@gmail.com',
    description='Mix up check tool',
    long_description=(open('README.md').read()),
    keywords='bioinformatics',
    url='https://github.com/AstraZeneca-NGS/ClearUp',
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
