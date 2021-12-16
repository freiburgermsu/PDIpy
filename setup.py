# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open('README.rst') as file:
    readme = file.read()

setup(
  name = 'PDIpy',      
  package_dir = {'core':'pdipy'},
  packages = find_packages(),
  package_data = {'pdipy': ['parameters/*.json']}, 
  version = '0.0.1',
  license = 'GNU',
  description = "A tool for calculating cellular oxidation and death via photodynamic inactivation.", 
  long_description = readme,
  long_description_content_type = "text/markdown",
  author = 'Andrew Freiburger',               
  author_email = 'andrewfreiburger@gmail.com',
  url = 'https://github.com/freiburgermsu/PDIpy',   
#   download_url = 'https://github.com/user/reponame/archive/v_01.tar.gz', 
  keywords = ['antibacterial', 'photodynamic', 'biophysics'],
  install_requires = ['matplotlib', 'tellurium', 'chemicals', 'scipy', 'pandas', 'sigfig']
)