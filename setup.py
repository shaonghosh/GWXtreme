import sys
from setuptools import setup

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup_requires = ['setuptools >= 30.3.0']
if {'build_sphinx'}.intersection(sys.argv):
    setup_requires.extend(['recommonmark', 'sphinx'])

setup(
     setup_requires=setup_requires,
     long_description=long_description,
     long_description_content_type='text/markdown'
)
