from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
  name = 'PDIpy',      
  packages = find_packages(where = "core"),
  version = '0.0.1',
  license = license,
  description = "A tool for calculating cellular oxidation and death via photodynamic inactivation.", 
  long_description = readme,
  long_description_content_type = "text/x-rst",
  author = 'Andrew Freiburger',               
  author_email = 'andrewfreiburger@gmail.com',
  url = 'https://github.com/freiburgermsu/PDIpy',   
  # download_url = 'https://github.com/user/reponame/archive/v_01.tar.gz', 
  keywords = ['antibacterial', 'photodynamic', 'biophysics'],   
  install_requires=[            
          'chemicals',
          'scipy',
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Developers and Users',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: GNU License',
    'Programming Language :: Python :: >=3',
  ],
)