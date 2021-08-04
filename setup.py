import setuptools
import sys

requirements = {"install": ["numpy"]}

if sys.version_info < (3, 7):
    requirements["install"].append("dataclasses")

version = open("xline/__init__.py").readline().split('"')[1]

setuptools.setup(
    name="xline",
    version=version,
    description="6D Tracking Code",
    author="Riccardo De Maria, Giovanni Iadarola",
    author_email="riccardo.de.maria@cern.ch",
    url="https://github.com/rdemaria/xline",
    packages=["xline", "xline.be_beamfields"],
    package_dir={"xline": "xline"},
    install_requires=requirements["install"],
)
