import setuptools

setuptools.setup(
        name='pysixtrack',
        version='0.0.1',
        description='6D Tracking Code',
        author='Riccardo De Maria',
        author_email='riccardo.de.maria@cern.ch',
        url='https://github.com/rdemaria/pysixtrack',
        packages=['pysixtrack'],
        package_dir={'pysixtrack': 'pysixtrack'},
        install_requires=['numpy']
)
