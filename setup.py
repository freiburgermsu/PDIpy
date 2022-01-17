# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open('README.rst') as file:
    readme = file.read()

setup(
  name = 'PDIpy',      
  package_dir = {'core':'pdipy'},
  packages = find_packages(),
  package_data = {'pdipy': ['parameters/*']}, 
  version = '0.0.2',
  license = 'MIT',
  description = "Predicted %-inactivation from a chemical kinetics model of cytoplasmic oxidation via photodynamic inactivation.", 
  long_description = readme,
  author = 'Andrew Freiburger',               
  author_email = 'andrewfreiburger@gmail.com',
  url = 'https://github.com/freiburgermsu/PDIpy',   
  keywords = ['antibacterial', 'photodynamic', 'biophysics', 'computational biology', 'medicine', 'PDI', 'antibiotics'],
  install_requires = ['matplotlib', 'tellurium', 'chemicals', 'scipy', 'pandas', 'sigfig', 'hillfit']
)