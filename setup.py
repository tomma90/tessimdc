import os
from setuptools import setup, find_packages

# Utility function to read the README file.
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = 'tessimdc',
    version = '0.0.1',
    author = 'Tommaso Ghigna',
    author_email = 't.ghigna@gmail.com',
    description = 'Simulation tool for DC biased TES detectors',
    license = 'GPLv3',
    keywords = 'TES detectors CMB cosmology',
    url = 'https://github.com/tomma90/tessimdc',
    packages = find_packages(),
    include_package_data=True,
    long_description = read('README.rst'),
    test_suite = 'tessimdc.test'
)
