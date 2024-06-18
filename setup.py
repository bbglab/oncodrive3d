from os import path

from setuptools import find_packages, setup

from scripts import __version__

include_package_data=True

directory = path.dirname(path.abspath(__file__))
with open(path.join(directory, 'requirements.txt')) as f:
    required = f.read().splitlines()

with open(path.join(directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='Oncodrive3D',
    version=__version__,
    description='BBGLab tool',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(),
    install_requires=required,
    include_package_data=True,
    url="https://github.com/bbglab/clustering_3d",
    author="BBGLab",
    author_email="stefano.pellegrini@irbbarcelona.org",
    license="GNU Affero General Public License v3 or later (AGPLv3+)", # TODO change accordingly
    entry_points={
        'console_scripts': [
            'oncodrive3D = scripts.main:oncodrive3D',
            'oncodrive3d = scripts.main:oncodrive3D',
            'Oncodrive3d = scripts.main:oncodrive3D',
            'Oncodrive3D = scripts.main:oncodrive3D'
        ]
    }
)