import sys
from setuptools import setup

setup_requires = ['setuptools >= 30.3.0']
if {'build_sphinx'}.intersection(sys.argv):
    setup_requires.extend(['recommonmark', 'sphinx'])

setup(setup_requires=setup_requires)
