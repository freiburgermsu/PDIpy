# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open('README.rst', 'r', encoding='utf-8') as file:
    readme = file.read()

setup(
  name = 'PDIpy',      
  package_dir = {'pdi':'pdipy'},
  packages = find_packages(),
  package_data = {
          'pdipy': ['parameters/*'],
          'test': ['./*']
          }, 
  version = '0.0.2',
  license = 'MIT',
  description = "Simulate Photodynamic Inactivation (PDI) of a Cocci Bacterium from a kinetics model of membrane oxidation.", 
  long_description = readme,
  author = 'Andrew Freiburger',               
  author_email = 'andrewfreiburger@gmail.com',
  url = 'https://github.com/freiburgermsu/PDIpy',   
  keywords = [
          'antibacterial',
          'photodynamic', 
          'biophysics',
          'computational',
          'biology',
          'medicine', 
          'PDI', 
          'antibiotics'
          ],
  install_requires = [
          'matplotlib',
          'tellurium', 
          'scipy', 
          'pandas',
          'sigfig',
          'hillfit',
          'chemw',
          'numpy'
          ]
)