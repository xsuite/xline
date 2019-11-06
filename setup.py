import setuptools

setuptools.setup(
    name="pysixtrack",
    version="0.0.3",
    description="6D Tracking Code",
    author="Riccardo De Maria",
    author_email="riccardo.de.maria@cern.ch",
    url="https://github.com/rdemaria/pysixtrack",
    packages=["pysixtrack","pysixtrack.be_beamfields"],
    package_dir={"pysixtrack": "pysixtrack"},
    install_requires=["numpy"],
)
