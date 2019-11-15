import setuptools
import sys

requirements={
    "install": ['numpy']
}
 
if sys.version_info < (3,7):
    requirements['install'].append('dataclasses')
    
setuptools.setup(
    name="pysixtrack",
    version="0.0.4",
    description="6D Tracking Code",
    author="Riccardo De Maria",
    author_email="riccardo.de.maria@cern.ch",
    url="https://github.com/rdemaria/pysixtrack",
    packages=["pysixtrack", "pysixtrack.be_beamfields"],
    package_dir={"pysixtrack": "pysixtrack"},
    install_requires=requirements['install'],
)
